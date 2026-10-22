// Pulse production for the AMoRE WCMD FADC raw data (NKFADC500, 2 ns/sample).
// Fills two trees: `normal` (charge = sum of found pulses in the signal window,
// dark_q in the front) and `pedestal` (0.1 Hz random trigger: largest
// excursion, window charges around it, dark_q over the record).
// Algorithm, settings and output: README.md.
//
//   produce_fadc -c prod_fadc.yaml -o out/prod_FADC1_<run>.<sub>.root
//                /data/amore2/data/RAW/<run>/FADC1_<run>.root.<sub>
//
// Several inputs (same run and board, FADC<b>_<run>.root.<subrun>) are chained
// into one output; --split rolls it over past that many bytes.

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

// GetID() is the PID, GetBit() the DAQ per-channel discriminator.
// Built-in defaults (DLY 1100); prod_fadc.yaml overrides them.
constexpr int kBaseLoDef = 0; // baseline region [baseLo, baseHi), in front of the signal
constexpr int kBaseHiDef = 208;
constexpr int kIntLoDef = 232; // integration region [intLo, intHi)
constexpr int kIntHiDef = 416;

constexpr double kPulseSigma = 4.0; // pulse-finding threshold, in baseline RMS
constexpr int kPulsePre = 3;        // samples added before each excursion
constexpr int kPulsePost = 21;      // samples added after each excursion
constexpr double kClipSigma = 4.0;  // baseline outlier rejection

constexpr Long64_t kSplitBytes = 1000000000LL; // a sensible --split value, not a default
constexpr int kMaxCh = 48;
constexpr int kMaxWin = 8; // integration windows the pedestal tree can carry

// ── Configuration ──────────────────────────────────────────────────────────
// Built-in defaults; the YAML overrides them and the values used go into the output.
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

  int pidLo = 0, pidHi = 0; // 0,0 accepts whatever the board sends

  TString path = "(built-in defaults)";
};

// Keys left out keep the default.  FADC1 keeps PID 1..47, FADC2 the door PMT 49.
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
// Even and odd samples separately (fixed parity offset), with one-sided
// clipping of pulses; [skipLo, skipHi) is left out.
struct Baseline {
  double even = 0.0;
  double odd = 0.0;
  double rms = 0.0; // after parity subtraction
};

Baseline EstimateBaseline(const unsigned short * wave, int lo, int hi, double clipSigma, int iters,
                          int skipLo = 0, int skipHi = 0)
{
  Baseline bl;
  // clip on the parity-subtracted deviation
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
// Excursions above thr in [lo, hi), widened by (pre, post), merged, summed;
// clamped to [lo, hi).
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
// One pass over the raw files fills both trees.
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
  printf("dark     : pulses in the baseline region, normal and pedestal\n");
  printf("pedestal : baseline %s minus [peak%+d,peak%+d)  peak in [%d,%d)  windows",
         cfg.pedHi > 0 ? Form("[%d,%d)", cfg.pedLo, cfg.pedHi) : "the whole record", cfg.exclLo,
         cfg.exclHi, cfg.hmaxLo, cfg.hmaxHi);
  for (const auto & w : cfg.windows)
    printf(" q_%d_%d", w.first, w.second);
  printf("\n");
  if (cfg.pidHi > 0) printf("channels : PID %d..%d\n", cfg.pidLo, cfg.pidHi);

  // Shared quantities have the same name and type in both trees.

  // ── normal ──
  Int_t runNo = run, subrun = 0, nch = 0, nhit = 0;
  UInt_t trgNum = 0, evtNum = 0;
  ULong64_t trgTime = 0;
  Short_t pid[kMaxCh], hmaxx[kMaxCh], npul[kMaxCh];
  UShort_t tbit[kMaxCh];
  Float_t base_even[kMaxCh], base_odd[kMaxCh], rms[kMaxCh];
  Float_t charge[kMaxCh], qwin[kMaxCh], hmax[kMaxCh], dark_q[kMaxCh];

  // ── pedestal ──
  Int_t p_subrun = 0, p_nch = 0;
  UInt_t p_trgNum = 0, p_evtNum = 0;
  ULong64_t p_trgTime = 0;
  Short_t p_pid[kMaxCh], p_hmaxx[kMaxCh], p_gap2[kMaxCh], p_qppos[kMaxCh];
  Float_t p_hmax[kMaxCh], p_hmax2[kMaxCh], p_be[kMaxCh], p_bo[kMaxCh], p_rms[kMaxCh];
  Float_t p_dark[kMaxCh];
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
    tnorm->Branch("dark_q", dark_q,
                  "dark_q[nch]/F"); // ADC x samples, sum of pulses in the baseline region

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
    tped->Branch("dark_q", p_dark,
                 "dark_q[nch]/F"); // ADC x samples, as the normal tree's, whole region
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
    TNamed("dark", "sum of the pulses (as `pulse`) in the baseline region").Write();
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

        // dark: the pulses in the baseline region, which holds no light
        const Charge dark =
            SumPulses(wave, bl, cfg.normLo, cfg.normHi, cfg.nsigma * bl.rms, cfg.pre, cfg.post);

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
        dark_q[nch] = static_cast<Float_t>(dark.q);
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

        // level from the whole record, then again without the peak
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
        p_dark[p_nch] = static_cast<Float_t>(
            SumPulses(wave, bl, plo, phi, cfg.nsigma * bl.rms, cfg.pre, cfg.post).q);
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
