// Pulse production for the AMoRE WCMD FADC raw data.
//
// For every waveform: estimate the baseline from the pedestal region at the
// front of the readout window, subtract it, and integrate the whole
// baseline-subtracted signal region to get the total charge in ADC counts x
// samples (1 sample = 2 ns, NKFADC500).
//
// Usage — raw files in, one output file out (prod.sh decides both for a run):
//   produce_fadc -c prod_fadc.yaml -o out/prod_FADC1_000263.00001.root
//                /data/amore2/data/RAW/000263/FADC1_000263.root.00001
//
// Several inputs are chained into one output, in the order given; they must
// all be the same run and board.  Run, board and subrun are read from the DAQ's
// file names, FADC<b>_<run>.root.<subrun>.
//
// --split splits the output, rolling over to <output>_1.root, _2, ... each
// time a file passes that many bytes; kSplitBytes (1 GB) is a sensible value.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <getopt.h>
#include <string>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "EventInfo.hh"
#include "FChannel.hh"
#include "FChannelData.hh"
#include "common.hh"
#include "yaml-cpp/yaml.h"

// ── Channel identity ───────────────────────────────────────────────────────
// FChannel::GetID() in these files is the PID, not the hardware CID: the DAQ
// writes `chdata->Add(conf->PID(i), ...)` (CUPDAQ_v4.0
// DAQ/src/CupDAQManager_write_root.cc:94).  CID repeats per FADC module
// (fadct.yml: CID [1,2,3,4] on every module, PID [1,2,3,4], [5,6,7,8], ...)
// and never reaches the ROOT file, so no remapping is needed or possible here.
//
// FChannel::GetBit() is the per-channel trigger bit from the ADC header
// (CupDAQManager_write_root.cc:98), i.e. the DAQ's own discriminator decision.
//
// Run 263: FADC1 carries PID 1-48, FADC2 carries PID 49 alone.  PID 48 is not
// connected — pure noise, trigger bit never set — but it is kept in the output
// rather than silently dropped; cut on `tbit` or `pid` downstream.

// ── Analysis windows ───────────────────────────────────────────────────────
// Three bands, all measured on run 263.  Note the 24-sample guard band: the
// pedestal ends at 208 and the integration starts at 232, they do not meet.
//
//     0 ............ 208 ... 232 ............... 416 ........... 1008
//     |   pedestal    | guard |   integration     |   nothing
//
// PEDESTAL [0, 208).  The fraction of samples above 5 sigma is at the
// accidental floor (0.004-0.04%) throughout, in all 24 subruns, and the signal
// edge sits at sample 232 in every one of them -- the trigger position does
// not drift over the 24 hour run.  Shortening this region only costs
// resolution (pedestal width 37.7 at [0,208), 41.5 at [0,150), 49.4 at [0,60))
// without buying margin, and does not shift the mean.  The ~1% of triggered
// channels with a 5 sigma crossing before 208 are uniform in time, i.e.
// accidental dark hits, not early light.
//
// GUARD [208, 232).  Belongs to neither: it is the leading edge of the signal,
// but it carries almost nothing.  Splitting its content against a signal-free
// reference region (to remove the baseline offset that comes from selecting
// channels on their own charge) leaves, as real light:
//
//   pulse size     50-100   100-200   200-400    >400
//   real charge      0.30      0.31      0.95   <0.40   ADC*sample
//   as %% of total    0.4%%     0.25%%     0.36%%   0.04%%
//
// Integrating it would add under 0.4% of signal and 24 samples of noise
// (pedestal width 33.5 -> 35.7), so the integration starts at the measured
// edge instead, and the band absorbs any small drift of that edge.
//
// INTEGRATION [232, 416).  The light is not prompt: photoelectrons trickle in
// over ~180 samples (~360 ns) because of the low light level and reflections
// in the tank, so the window has to span the whole structure rather than sit
// on the peak.  Mean cumulative charge, normalised to its own row maximum, for
// triggered tank channels split by pulse size:
//
//   end sample     304    352    384    416    448    512    672
//   q =  50-100   0.802  0.970  0.981  0.969  0.950  0.920  0.931
//   q = 100-200   0.582  0.741  0.773  0.780  0.778  0.780  0.835
//   q = 200-400   0.596  0.906  0.994  1.000  0.981  0.940  0.924
//   q >     400   0.236  0.862  0.988  1.000  0.977  0.913  0.856
//
// Large pulses have proportionally longer tails, so a short window is not
// merely a scale factor: cutting at 352 costs 3% on small pulses but 14% on
// the largest, a real non-linearity.  Past ~416 the AC-coupling undershoot
// eats the integral back (q > 400 loses 14% by sample 672), so a longer window
// both adds noise and removes signal.  Between 384 and 416 the charge agrees
// to within 1.2% across the dynamic range; 416 is where the two large-pulse
// classes peak.
//
// The remaining cost is resolution: the pedestal width of the charge grows
// roughly linearly with window length (~0.18 ADC per sample: 5.7 ADC*sample
// over 20 samples, 33.5 over 184, 35.7 over 208), because the baseline wanders
// faster than the per-sample white noise of ~1.04 ADC would suggest.  Taking
// the baseline from both sides of the signal instead of only the front buys
// just 14%, so window length is the only real lever -- hence these are
// parameters, not constants, so they can be scanned.

constexpr int kBaseLoDef = 0; // pedestal region [baseLo, baseHi)
constexpr int kBaseHiDef = 208;
constexpr int kIntLoDef = 232; // integration region [intLo, intHi)
constexpr int kIntHiDef = 416;

// ── Trigger type ───────────────────────────────────────────────────────────
// trgType 1 is the pedestal (random) trigger and is skipped by default.
// Measured on run 263: consecutive type-1 triggers are exactly 10 s apart
// (min = median = max = 1e10 ns) on both boards independently, i.e. a 0.1 Hz
// clock.  They carry no light — mean charge 0.37 vs 11.14 ADC*sample on FADC1,
// signal-window occupancy 0.004% vs 0.32%, trigger bit set on 0.2% of channels
// vs 11.4% — and on FADC2, which self-triggers, no channel has the trigger bit
// set at all.  They are 0.5% of the events.
//
// Only types 0 and 1 occur in the whole run: FADC1 1669186 / 8640, FADC2
// 2531109 / 8639.  8640 = 24 h / 10 s exactly.  Neither tree carries trgType:
// which tree an event went to says which it was.  kPedTrigType is in common.hh.

// ── Charge ─────────────────────────────────────────────────────────────────
// The charge is the sum of the found pulses, not the integral of the whole
// window.  The Cherenkov signal is a sparse train of narrow pulses, so most of
// a fixed window is empty and only contributes noise: measured on a
// all 8640 pedestal-trigger events of run 263 -- which contain no light at all
// -- the response is
//
//                      per channel          summed over the 47 tank channels
//   fixed window   mean 0.39  RMS 31.73     mean 18.27  RMS 268.27
//   sum of pulses  mean 0.30  RMS  3.66     mean 14.28  RMS  30.39
//
// 97.1% of channels give exactly zero.  At event level, where it matters for
// the muon light yield, the noise floor drops from 268 to 30 ADC*sample --
// 8.8x -- against a typical event charge of ~500.  The small residual bias
// (14.3 per event, from the 2.9% of channels that catch a noise excursion) is
// constant and can be subtracted.  The single-pulse template (aligned on its peak,
// single-pulse channels only) is narrow — offsets -3..+5 already hold 86% of
// the pulse, the rest is a low tail — so widening each excursion by kPulsePre
// before and kPulsePost after recovers essentially all of it.  Charge retained
// relative to the fixed window, by number of pulses found:
//
//   npul     1      2      3      4      5      6      7      8
//   ratio  0.999  0.978  0.983  0.995  0.999  1.005  1.000  0.998
//
// i.e. within 2.2% across the whole dynamic range, so the gain in resolution
// costs almost no linearity.  Below 4 sigma the threshold starts catching
// noise (pedestal RMS 5.6 at 3 sigma, 9.5 at 2.5 sigma); above it, small
// pulses are missed.  post=21 is where the ratio reaches 1 (0.92 at post=5,
// 0.97 at 13, 1.02 at 33 where it starts picking up baseline).  `pre` barely
// matters: the template holds only 1.3% of the charge before peak-3.
//
// The fixed-window integral is kept alongside as `qwin` for cross-checks.

constexpr double kPulseSigma = 4.0; // pulse-finding threshold, in baseline RMS
constexpr int kPulsePre = 3;        // samples added before each excursion
constexpr int kPulsePost = 21;      // samples added after each excursion

// ── PMT hit ────────────────────────────────────────────────────────────────
// A hit is the DAQ trigger bit AND a charge above kHitCharge.  The trigger bit
// is a clean height discriminator: P(tbit) is 0.004 below hmax 14 ADC, 0.66 at
// 14-16, 0.99 at 16-18 and exactly 1.0000 above 18, identically on all 47 tank
// channels, and no channel with hmax > 25 in the window ever has tbit = 0.
// (The 3.7% of tbit = 1 channels with no pulse in the window are dark pulses
// elsewhere in the 1008-sample record, which the discriminator also watches.)
//
// The charge cut sits above the noise peak at ~8 and the valley at ~19 in the
// npul > 0 spectrum.  On top of tbit it removes only ~0.1% more, so it is
// there to make the definition independent of the trigger bit rather than to
// clean the sample.
//
// Note for gain work: tbit already truncates the low-charge side of the SPE
// peak - its efficiency is 0.21 at charge 30-40, 0.64 at 40-50, 0.87 at 50-58 -
// so the SPE mode and mean come out 5% and 9% high compared with `npul > 0 &&
// charge > 25`.  Use that selection instead when calibrating gain.

constexpr double kHitCharge = 25.0; // ADC x samples
constexpr double kClipSigma = 4.0;  // pedestal outlier rejection

// Splitting the output is off by default: `maxFileBytes <= 0` writes one file
// however large it grows.  Passing a size rolls the output over the way ROOT's
// own TTree::SetMaxTreeSize would, except that each file keeps its own copy of
// the metadata below rather than only the last one.  This is a suggested value
// to pass, not a default.

constexpr Long64_t kSplitBytes = 1000000000LL; // 1 GB
constexpr int kMaxCh = 48;
constexpr int kMaxWin = 8; // integration windows the pedestal tree can carry

// ── Configuration ──────────────────────────────────────────────────────────
// Everything above is a default.  What a run actually used comes from a YAML
// file, because these numbers are still being studied and a production is only
// reproducible if the settings travel with it rather than living in a rebuild.
struct Config {
  double clipSigma = kClipSigma;
  int iters = 4;

  int normLo = kBaseLoDef, normHi = kBaseHiDef; // normal: region the level comes from
  int pedLo = 0, pedHi = 0;                     // pedestal: same, 0 for the whole record
  int exclLo = -15, exclHi = 50;                // pedestal: cut out around the peak

  int intLo = kIntLoDef, intHi = kIntHiDef;
  double nsigma = kPulseSigma;
  int pre = kPulsePre, post = kPulsePost;

  std::vector<std::pair<int, int>> windows{{2, 6}, {4, 8}, {6, 16}, {8, 24}, {10, 45}};
  int hmaxLo = 15, hmaxHi = 957;
  int pedAway = 80, secondAway = 20;
  int seed = 0;

  // FADC1 writes 48 channels but only 1..47 carry a PMT; the 48th PMT is at
  // the door and reads out on FADC2 as PID 49, on its own trigger path.
  int pidLo = 0, pidHi = 0; // 0,0 accepts whatever the board sends

  TString path = "(built-in defaults)";
};

// Read what the file gives and leave the rest alone, so a partial config is a
// set of overrides rather than a trap.
bool LoadConfig(const char * path, Config & c, int board)
{
  if (board == 1) {
    c.pidLo = 1;
    c.pidHi = 47;
  }
  else if (board == 2) {
    c.pidLo = 49;
    c.pidHi = 49;
  }
  if (!path || !path[0]) return true;
  YAML::Node n;
  if (!LoadYaml(path, n)) return false;
  const auto pair = ReadPair;

  if (const YAML::Node b = n["baseline"]) {
    if (b["clip_sigma"]) c.clipSigma = b["clip_sigma"].as<double>();
    if (b["iters"]) c.iters = b["iters"].as<int>();
    pair(b["window_normal"], c.normLo, c.normHi);
    pair(b["window_pedestal"], c.pedLo, c.pedHi);
    pair(b["exclude_pedestal"], c.exclLo, c.exclHi);
  }
  if (const YAML::Node m = n["normal"]) {
    pair(m["integration"], c.intLo, c.intHi);
    if (const YAML::Node p = m["pulse"]) {
      if (p["nsigma"]) c.nsigma = p["nsigma"].as<double>();
      if (p["pre"]) c.pre = p["pre"].as<int>();
      if (p["post"]) c.post = p["post"].as<int>();
    }
  }
  if (const YAML::Node d = n["pedestal"]) {
    if (const YAML::Node w = d["windows"]) {
      c.windows.clear();
      for (const auto & e : w)
        c.windows.emplace_back(e[0].as<int>(), e[1].as<int>());
    }
    pair(d["hmax_range"], c.hmaxLo, c.hmaxHi);
    if (d["ped_away"]) c.pedAway = d["ped_away"].as<int>();
    if (d["second_away"]) c.secondAway = d["second_away"].as<int>();
    if (d["seed"]) c.seed = d["seed"].as<int>();
  }
  if (const YAML::Node k = n["channels"]) {
    if (const YAML::Node b = k[board]) pair(b, c.pidLo, c.pidHi);
  }
  c.path = path;
  return true;
}

// ── Baseline ───────────────────────────────────────────────────────────────
// Every channel shows a fixed even/odd sample offset (up to +-1.3 ADC on the
// door PMT, PID-dependent on the tank channels), so the baseline is estimated
// separately for even and odd samples and subtracted per parity.  A few
// iterations of one-sided clipping keep a stray dark pulse inside the region
// from pulling the estimate up.
//
// [skipLo, skipHi) is left out, which is how the pedestal trigger can take its
// baseline from the whole record apart from the pulse it just found.  A
// pedestal-trigger record has no signal region, so using all of it rather than
// only the front measurably tightens the charge: on run 263 the pedestal width
// of a 13-sample window is 3.77 against 4.13, 9.5% better, because the level
// wanders over the record faster than white noise would.
struct Baseline {
  double even = 0.0;
  double odd = 0.0;
  double rms = 0.0; // after parity subtraction
};

Baseline EstimateBaseline(const unsigned short * wave, int lo, int hi, double clipSigma, int iters,
                          int skipLo = 0, int skipHi = 0)
{
  Baseline bl;
  // The rejection is on the parity-subtracted deviation, not on the raw sample
  // against the parity-averaged level: a channel whose two parities sit up to
  // 3 ADC apart would otherwise have one of them cut far harder than the other.
  double selEven = 0.0, selOdd = 0.0, selCut = 1e9;
  long selN = -1;

  for (int iter = 0; iter < iters; iter++) {
    auto kept = [&](int s) {
      if (s >= skipLo && s < skipHi) return false;
      if (selN < 0) return true; // nothing selected yet: take everything
      return wave[s] - (s % 2 == 0 ? selEven : selOdd) <= selCut;
    };

    double se = 0.0, so = 0.0;
    long ne = 0, no = 0;
    for (int s = lo; s < hi; s++) {
      if (!kept(s)) continue;
      if (s % 2 == 0) {
        se += wave[s];
        ne++;
      }
      else {
        so += wave[s];
        no++;
      }
    }
    if (ne == 0 || no == 0) break;
    bl.even = se / ne;
    bl.odd = so / no;

    double s2 = 0.0;
    long n = 0;
    for (int s = lo; s < hi; s++) {
      if (!kept(s)) continue;
      const double d = wave[s] - (s % 2 == 0 ? bl.even : bl.odd);
      s2 += d * d;
      n++;
    }
    if (n == 0) break;
    bl.rms = std::sqrt(s2 / n);

    // how many survive the new cut, taken with the level just measured
    const double cut = clipSigma * bl.rms;
    long nNew = 0;
    for (int s = lo; s < hi; s++) {
      if (s >= skipLo && s < skipHi) continue;
      if (wave[s] - (s % 2 == 0 ? bl.even : bl.odd) <= cut) nNew++;
    }
    if (nNew == selN) break; // converged
    selEven = bl.even;
    selOdd = bl.odd;
    selCut = cut;
    selN = nNew;
  }

  return bl;
}

inline double Sub(const unsigned short * wave, const Baseline & bl, int s)
{
  return wave[s] - (s % 2 == 0 ? bl.even : bl.odd);
}

// ── Pulse finding ──────────────────────────────────────────────────────────
// Find the excursions above threshold in [lo, hi), widen each by (pre, post),
// merge those that then overlap, and sum the baseline-subtracted samples they
// cover.  Windows are clamped to [lo, hi), which loses nothing: the signal is
// already back at the noise floor before hi.
struct Charge {
  double q = 0.0;
  int npulse = 0;
};

Charge SumPulses(const unsigned short * wave, const Baseline & bl, int lo, int hi, double thr,
                 int pre, int post)
{
  Charge c;
  int wStart = 0, wEnd = -1; // current merged window, wEnd < wStart means none open
  bool above = false;
  int aStart = 0;

  auto flush = [&]() {
    if (wEnd < wStart) return;
    for (int s = std::max(wStart, lo); s <= wEnd && s < hi; s++)
      c.q += Sub(wave, bl, s);
  };

  for (int s = lo; s <= hi; s++) {
    const bool over = (s < hi) && (Sub(wave, bl, s) > thr);
    if (over && !above) {
      above = true;
      aStart = s;
    }
    else if (!over && above) {
      above = false;
      c.npulse++;
      const int nStart = aStart - pre, nEnd = s + post;
      if (wEnd >= wStart && nStart <= wEnd + 1) { // overlaps the open window
        wEnd = std::max(wEnd, nEnd);
      }
      else {
        flush();
        wStart = nStart;
        wEnd = nEnd;
      }
    }
  }
  flush();

  return c;
}

// ── Producer ───────────────────────────────────────────────────────────────
// One pass over the raw file feeds both trees: `normal` for the normal
// triggers and `pedestal` for the 0.1 Hz pedestal ones.  They want different
// things -- the normal-trigger charge comes from the pulses found inside a fixed
// signal window, the pedestal charge from a set of windows around whatever the
// largest excursion in the record happens to be -- but they read the same
// waveform, so splitting them into two productions would read a terabyte
// twice.
int produce(const std::vector<std::string> & inputs, const char * output, Long64_t maxEvents,
            Long64_t maxFileBytes, const char * cfgPath)
{
  int run = -1, board = -1;
  std::vector<int> subrunOfTree; // which subrun each chained file came from
  for (const std::string & f : inputs) {
    TString stream;
    int r = -1, sub = -1;
    if (!ParseRawName(f.c_str(), stream, r, sub) || !stream.BeginsWith("FADC")) {
      printf("ERROR: '%s' is not named FADC<b>_<run>.root.<subrun>\n", f.c_str());
      return 1;
    }
    const int b = TString(stream(4, stream.Length())).Atoi();
    if (run < 0) {
      run = r;
      board = b;
    }
    else if (r != run || b != board) {
      printf("ERROR: %s is not run %d FADC%d like the first input\n", f.c_str(), run, board);
      return 1;
    }
    if (gSystem->AccessPathName(f.c_str())) {
      printf("ERROR: no input file %s\n", f.c_str());
      return 1;
    }
    subrunOfTree.push_back(sub);
  }

  Config cfg;
  if (!LoadConfig(cfgPath, cfg, board)) return 1;

  if (cfg.normHi <= cfg.normLo || cfg.intHi <= cfg.intLo) {
    printf("ERROR: empty window: baseline [%d,%d)  integration [%d,%d)\n", cfg.normLo, cfg.normHi,
           cfg.intLo, cfg.intHi);
    return 1;
  }
  if (cfg.normHi > cfg.intLo) printf("WARNING: baseline region overlaps the integration region\n");
  if ((cfg.intHi - cfg.intLo) % 2 != 0)
    printf("WARNING: integration window has an odd number of samples\n");
  if (cfg.windows.empty() || (int)cfg.windows.size() > kMaxWin) {
    printf("ERROR: %d pedestal windows, must be 1..%d\n", (int)cfg.windows.size(), kMaxWin);
    return 1;
  }
  if (board < 1) {
    printf("ERROR: board must be >= 1\n");
    return 1;
  }

  int maxHead = 0, maxTail = 0;
  for (const auto & w : cfg.windows) {
    maxHead = std::max(maxHead, w.first);
    maxTail = std::max(maxTail, w.second);
  }
  if (cfg.hmaxLo < maxHead || cfg.hmaxLo < -cfg.exclLo) {
    printf("ERROR: hmax_range starts at %d, too early for head %d / exclude %d\n", cfg.hmaxLo,
           maxHead, cfg.exclLo);
    return 1;
  }

  const int from = subrunOfTree.front(), to = subrunOfTree.back();
  const TString subLabel =
      from == to ? TString::Format("%d", from) : TString::Format("%d-%d", from, to);
  const TString inLabel = inputs.size() == 1
                              ? TString(inputs[0])
                              : TString::Format("%s .. %s (%zu files)", inputs.front().c_str(),
                                                inputs.back().c_str(), inputs.size());
  TString outBase = output;
  if (outBase.EndsWith(".root")) outBase.Resize(outBase.Length() - 5);

  auto * chain = new TChain("AbsEvent");
  for (const std::string & f : inputs)
    chain->Add(f.c_str());

  EventInfo * info = nullptr;
  FChannelData * data = nullptr;
  chain->SetBranchAddress("EventInfo", &info);
  chain->SetBranchAddress("FChannelData", &data);

  Long64_t nentries = chain->GetEntries();
  if (maxEvents > 0 && maxEvents < nentries) nentries = maxEvents;
  printf("run %d  subrun %s  board FADC%d\n", run, subLabel.Data(), board);
  printf("config : %s\n", cfg.path.Data());
  printf("input  : %s  (%lld events in %d files)\n", inLabel.Data(), nentries,
         chain->GetListOfFiles()->GetEntries());
  printf("output : %s  %s\n", output,
         maxFileBytes > 0 ? Form("(new file past %.2f GB)", maxFileBytes / 1e9)
                          : "(single file, no split)");
  printf("normal   : baseline [%d,%d)  integration [%d,%d)  pulses %.1f sigma (-%d,+%d)\n",
         cfg.normLo, cfg.normHi, cfg.intLo, cfg.intHi, cfg.nsigma, cfg.pre, cfg.post);
  printf("pedestal : baseline %s minus [peak%+d,peak%+d)  peak in [%d,%d)  windows",
         cfg.pedHi > 0 ? Form("[%d,%d)", cfg.pedLo, cfg.pedHi) : "the whole record", cfg.exclLo,
         cfg.exclHi, cfg.hmaxLo, cfg.hmaxHi);
  for (const auto & w : cfg.windows)
    printf(" q_%d_%d", w.first, w.second);
  printf("\n");
  if (cfg.pidHi > 0) printf("channels : PID %d..%d\n", cfg.pidLo, cfg.pidHi);

  // The two trees share every quantity they both have, under the same name and
  // the same type, so code that reads one reads the other.
  //
  //   per event    subrun evtNum trgNum trgTime nch      (+ nhit on normal)
  //   per channel  pid base_even base_odd rms hmax hmaxx
  //
  // What is left is what only one of them measures.

  // ── normal ──
  Int_t runNo = run, subrun = 0, nch = 0, nhit = 0;
  UInt_t trgNum = 0, evtNum = 0;
  ULong64_t trgTime = 0;
  Short_t pid[kMaxCh], hmaxx[kMaxCh], npul[kMaxCh];
  UShort_t tbit[kMaxCh];
  Float_t base_even[kMaxCh], base_odd[kMaxCh], rms[kMaxCh];
  Float_t charge[kMaxCh], qwin[kMaxCh], hmax[kMaxCh];

  // ── pedestal ──
  Int_t p_subrun = 0, p_nch = 0;
  UInt_t p_trgNum = 0, p_evtNum = 0;
  ULong64_t p_trgTime = 0;
  Short_t p_pid[kMaxCh], p_hmaxx[kMaxCh], p_gap2[kMaxCh], p_qppos[kMaxCh];
  Float_t p_hmax[kMaxCh], p_hmax2[kMaxCh], p_be[kMaxCh], p_bo[kMaxCh], p_rms[kMaxCh];
  Float_t p_q[kMaxWin][kMaxCh], p_qp[kMaxWin][kMaxCh];
  const int nwin = cfg.windows.size();

  TFile * fout = nullptr;
  TTree * tnorm = nullptr;
  TTree * tped = nullptr;
  int seq = 0;
  Long64_t nnorm = 0, nped = 0;

  auto openOut = [&]() {
    const TString name =
        seq == 0 ? TString(output) : TString::Format("%s_%d.root", outBase.Data(), seq);
    fout = TFile::Open(name, "recreate", "", 1);
    if (!fout || fout->IsZombie()) {
      printf("ERROR: cannot write '%s'\n", name.Data());
      fout = nullptr;
      return false;
    }

    tnorm = new TTree("normal", "WCMD normal-trigger production");
    tnorm->Branch("run", &runNo, "run/I");
    tnorm->Branch("nch", &nch, "nch/I");
    tnorm->Branch("subrun", &subrun, "subrun/I");
    tnorm->Branch("evtNum", &evtNum, "evtNum/i");
    tnorm->Branch("trgNum", &trgNum, "trgNum/i");
    tnorm->Branch("trgTime", &trgTime, "trgTime/l");           // ns since the run started
    tnorm->Branch("nhit", &nhit, "nhit/I");                    // channels with the trigger bit set
    tnorm->Branch("pid", pid, "pid[nch]/S");                   // channel id, as the DAQ writes it
    tnorm->Branch("base_even", base_even, "base_even[nch]/F"); // raw ADC
    tnorm->Branch("base_odd", base_odd, "base_odd[nch]/F");
    tnorm->Branch("rms", rms, "rms[nch]/F");          // ADC, after parity subtraction
    tnorm->Branch("hmax", hmax, "hmax[nch]/F");       // ADC over baseline, in the window
    tnorm->Branch("hmaxx", hmaxx, "hmaxx[nch]/S");    // where it is
    tnorm->Branch("tbit", tbit, "tbit[nch]/s");       // DAQ per-channel discriminator
    tnorm->Branch("charge", charge, "charge[nch]/F"); // ADC x samples, sum of found pulses
    tnorm->Branch("qwin", qwin, "qwin[nch]/F");       // ADC x samples, plain window integral
    tnorm->Branch("npul", npul, "npul[nch]/S");       // pulses found in the window

    tped = new TTree("pedestal", "WCMD pedestal-trigger production");
    tped->Branch("nch", &p_nch, "nch/I");
    tped->Branch("subrun", &p_subrun, "subrun/I");
    tped->Branch("evtNum", &p_evtNum, "evtNum/i");
    tped->Branch("trgNum", &p_trgNum, "trgNum/i");
    tped->Branch("trgTime", &p_trgTime, "trgTime/l");    // ns since the run started
    tped->Branch("pid", p_pid, "pid[nch]/S");            // channel id, as the DAQ writes it
    tped->Branch("base_even", p_be, "base_even[nch]/F"); // raw ADC
    tped->Branch("base_odd", p_bo, "base_odd[nch]/F");
    tped->Branch("rms", p_rms, "rms[nch]/F");         // ADC, after parity subtraction
    tped->Branch("hmax", p_hmax, "hmax[nch]/F");      // ADC over baseline, anywhere in range
    tped->Branch("hmaxx", p_hmaxx, "hmaxx[nch]/S");   // where it is
    tped->Branch("hmax2", p_hmax2, "hmax2[nch]/F");   // next largest, well away from the peak
    tped->Branch("gap2", p_gap2, "gap2[nch]/S");      // signed distance to it
    tped->Branch("qp_pos", p_qppos, "qp_pos[nch]/S"); // where the pulse-free gate was put
    for (int w = 0; w < nwin; w++) {
      const TString q = TString::Format("q_%d_%d", cfg.windows[w].first, cfg.windows[w].second);
      tped->Branch(q, p_q[w], TString::Format("%s[nch]/F", q.Data()));
      const TString qp = "qp_" + q(2, q.Length());
      tped->Branch(qp, p_qp[w], TString::Format("%s[nch]/F", qp.Data()));
    }

    TNamed("run", Form("%d", run)).Write();
    TNamed("board", Form("%d", board)).Write();
    TNamed("subrun", subLabel.Data()).Write();
    TNamed("input", inLabel.Data()).Write();
    TNamed("config", cfg.path.Data()).Write();
    TNamed("window_normal", Form("[%d,%d)", cfg.normLo, cfg.normHi)).Write();
    TNamed("window_pedestal",
           cfg.pedHi > 0 ? Form("[%d,%d)", cfg.pedLo, cfg.pedHi) : "whole record")
        .Write();
    TNamed("exclude_pedestal", Form("[%d,%d)", cfg.exclLo, cfg.exclHi)).Write();
    TNamed("clip_sigma", Form("%.2f", cfg.clipSigma)).Write();
    TNamed("iters", Form("%d", cfg.iters)).Write();
    TNamed("integration", Form("[%d,%d)", cfg.intLo, cfg.intHi)).Write();
    TNamed("pulse", Form("%.2f sigma, -%d +%d", cfg.nsigma, cfg.pre, cfg.post)).Write();
    TNamed("hmax_range", Form("[%d,%d)", cfg.hmaxLo, cfg.hmaxHi)).Write();
    TNamed("ped_away", Form("%d", cfg.pedAway)).Write();
    TNamed("second_away", Form("%d", cfg.secondAway)).Write();
    TNamed("seed", Form("%d", cfg.seed)).Write();
    TNamed("channels", Form("%d..%d", cfg.pidLo, cfg.pidHi)).Write();
    TNamed("maxFileBytes", Form("%lld", maxFileBytes)).Write();
    return true;
  };

  auto closeOut = [&](bool more) {
    if (!fout) return;
    fout->cd();
    nnorm += tnorm->GetEntries();
    nped += tped->GetEntries();
    printf("  wrote %s : %lld normal, %lld pedestal, %.2f GB\n", fout->GetName(),
           tnorm->GetEntries(), tped->GetEntries(), fout->GetBytesWritten() / 1e9);
    tnorm->Write();
    tped->Write();
    fout->Close();
    fout = nullptr;
    tnorm = tped = nullptr;
    if (more) seq++;
  };

  if (!openOut()) return 1;

  TRandom3 rng(cfg.seed + 1000 * run + from);
  TStopwatch watch;
  watch.Start();
  Long64_t nshort = 0;

  for (Long64_t i = 0; i < nentries; i++) {
    chain->GetEntry(i);
    const bool isPed = info->GetTriggerType() == kPedTrigType;

    if (!isPed) {
      subrun = subrunOfTree[chain->GetTreeNumber()];
      trgNum = info->GetTriggerNumber();
      evtNum = info->GetEventNumber();
      trgTime = info->GetTriggerTime();
      nch = 0;
      nhit = 0;

      for (int j = 0; j < data->GetN(); j++) {
        const FChannel * ch = data->Get(j);
        if (nch >= kMaxCh) {
          printf("WARNING: event %lld has more than %d channels\n", i, kMaxCh);
          break;
        }
        if (cfg.pidHi > 0 && (ch->GetID() < cfg.pidLo || ch->GetID() > cfg.pidHi)) continue;
        const unsigned short * wave = ch->GetWaveform();
        if (ch->GetNdp() < cfg.intHi) {
          nshort++;
          continue;
        }

        const Baseline bl =
            EstimateBaseline(wave, cfg.normLo, cfg.normHi, cfg.clipSigma, cfg.iters);

        double qfix = 0.0, amax = -1e9;
        int iamax = cfg.intLo;
        for (int s = cfg.intLo; s < cfg.intHi; s++) {
          const double d = Sub(wave, bl, s);
          qfix += d;
          if (d > amax) {
            amax = d;
            iamax = s;
          }
        }
        const Charge pulses =
            SumPulses(wave, bl, cfg.intLo, cfg.intHi, cfg.nsigma * bl.rms, cfg.pre, cfg.post);

        pid[nch] = static_cast<Short_t>(ch->GetID());
        tbit[nch] = ch->GetBit();
        base_even[nch] = static_cast<Float_t>(bl.even);
        base_odd[nch] = static_cast<Float_t>(bl.odd);
        rms[nch] = static_cast<Float_t>(bl.rms);
        charge[nch] = static_cast<Float_t>(pulses.q);
        qwin[nch] = static_cast<Float_t>(qfix);
        npul[nch] = static_cast<Short_t>(pulses.npulse);
        hmax[nch] = static_cast<Float_t>(amax);
        hmaxx[nch] = static_cast<Short_t>(iamax);
        nch++;

        if (ch->GetBit() != 0) nhit++;
      }
      tnorm->Fill();
    }
    else {
      p_subrun = subrunOfTree[chain->GetTreeNumber()];
      p_evtNum = info->GetEventNumber();
      p_trgNum = info->GetTriggerNumber();
      p_trgTime = info->GetTriggerTime();
      p_nch = 0;

      for (int j = 0; j < data->GetN(); j++) {
        const FChannel * ch = data->Get(j);
        if (p_nch >= kMaxCh) break;
        if (cfg.pidHi > 0 && (ch->GetID() < cfg.pidLo || ch->GetID() > cfg.pidHi)) continue;
        const unsigned short * wave = ch->GetWaveform();
        const int ndp = ch->GetNdp();
        if (ndp < cfg.hmaxHi + maxTail + 1) {
          nshort++;
          continue;
        }

        // the whole record is pedestal here, so take the level from all of it,
        // then again with the pulse it found cut out
        const int plo = cfg.pedLo;
        const int phi = cfg.pedHi > 0 ? std::min(cfg.pedHi, ndp) : ndp;
        Baseline bl = EstimateBaseline(wave, plo, phi, cfg.clipSigma, cfg.iters);
        int hx = cfg.hmaxLo;
        double best = -1e9;
        for (int s = cfg.hmaxLo; s < cfg.hmaxHi; s++) {
          const double d = Sub(wave, bl, s);
          if (d > best) {
            best = d;
            hx = s;
          }
        }
        bl = EstimateBaseline(wave, plo, phi, cfg.clipSigma, cfg.iters, hx + cfg.exclLo,
                              hx + cfg.exclHi);
        best = -1e9;
        for (int s = cfg.hmaxLo; s < cfg.hmaxHi; s++) {
          const double d = Sub(wave, bl, s);
          if (d > best) {
            best = d;
            hx = s;
          }
        }

        // the same integrator somewhere with no pulse in it
        int pp = hx;
        for (int t = 0; t < 8; t++) {
          pp = cfg.hmaxLo + static_cast<int>(rng.Uniform(cfg.hmaxHi - cfg.hmaxLo));
          if (std::abs(pp - hx) >= cfg.pedAway) break;
        }

        // the next largest thing in the record, whatever its size
        double second = -1e9;
        int i2 = hx;
        for (int s = cfg.hmaxLo; s < cfg.hmaxHi; s++) {
          if (std::abs(s - hx) < cfg.secondAway) continue;
          const double d = Sub(wave, bl, s);
          if (d > second) {
            second = d;
            i2 = s;
          }
        }

        p_pid[p_nch] = static_cast<Short_t>(ch->GetID());
        p_hmaxx[p_nch] = static_cast<Short_t>(hx);
        p_gap2[p_nch] = static_cast<Short_t>(i2 - hx);
        p_qppos[p_nch] = static_cast<Short_t>(pp);
        p_hmax[p_nch] = static_cast<Float_t>(best);
        p_hmax2[p_nch] = static_cast<Float_t>(second);
        p_be[p_nch] = static_cast<Float_t>(bl.even);
        p_bo[p_nch] = static_cast<Float_t>(bl.odd);
        p_rms[p_nch] = static_cast<Float_t>(bl.rms);
        for (int w = 0; w < nwin; w++) {
          const int hd = cfg.windows[w].first, tl = cfg.windows[w].second;
          double q = 0.0, qp = 0.0;
          for (int s = hx - hd; s <= hx + tl; s++)
            q += Sub(wave, bl, s);
          for (int s = pp - hd; s <= pp + tl; s++)
            qp += Sub(wave, bl, s);
          p_q[w][p_nch] = static_cast<Float_t>(q);
          p_qp[w][p_nch] = static_cast<Float_t>(qp);
        }
        p_nch++;
      }
      tped->Fill();
    }

    if (maxFileBytes > 0 && fout->GetBytesWritten() > maxFileBytes) {
      closeOut(true);
      if (!openOut()) return 1;
    }

    if (i > 0 && i % 10000 == 0) {
      printf("  %lld / %lld  (%.1f s)\n", i, nentries, watch.RealTime());
      watch.Continue();
    }
  }

  closeOut(false);

  if (nshort > 0) printf("WARNING: %lld waveforms were too short and were skipped\n", nshort);
  printf("done: %lld events read -> %lld normal + %lld pedestal, %d file%s, in %.1f s\n", nentries,
         nnorm, nped, seq + 1, seq ? "s" : "", watch.RealTime());
  return 0;
}

// ── Command line ───────────────────────────────────────────────────────────
static void Usage(const char * prog)
{
  printf("usage: %s -o OUTPUT [-c CONFIG] [-n MAXEVENTS] [--split BYTES] INPUT...\n"
         "  INPUT...            raw FADC<b>_<run>.root.<subrun> files, chained in this order\n"
         "  -o, --out FILE      output file (required)\n"
         "  -c, --config FILE   YAML settings (default: built-in)\n"
         "  -n, --max-events N  stop after N events\n"
         "      --split BYTES   new output file past BYTES (%lld = 1 GB is sensible)\n",
         prog, kSplitBytes);
}

int main(int argc, char ** argv)
{
  long long maxEvents = -1, split = 0;
  const char *cfg = "", *out = nullptr;

  const option opts[] = {{"out", required_argument, nullptr, 'o'},
                         {"config", required_argument, nullptr, 'c'},
                         {"max-events", required_argument, nullptr, 'n'},
                         {"split", required_argument, nullptr, 'S'},
                         {"help", no_argument, nullptr, 'h'},
                         {nullptr, 0, nullptr, 0}};
  int c;
  while ((c = getopt_long(argc, argv, "o:c:n:h", opts, nullptr)) != -1) {
    bool ok = true;
    switch (c) {
      case 'o': out = optarg; break;
      case 'c': cfg = optarg; break;
      case 'n': ok = ParseInt(optarg, maxEvents); break;
      case 'S': ok = ParseInt(optarg, split); break;
      case 'h': Usage(argv[0]); return 0;
      default: Usage(argv[0]); return 2;
    }
    if (!ok) {
      printf("ERROR: bad value '%s' for %s\n", optarg, c == 'S' ? "--split" : "-n");
      return 2;
    }
  }
  if (!out || optind == argc) {
    Usage(argv[0]);
    return 2;
  }
  const std::vector<std::string> inputs(argv + optind, argv + argc);
  return produce(inputs, out, maxEvents, split, cfg);
}
