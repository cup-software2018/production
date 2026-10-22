// Pulse production for the AMoRE PSMD (NKIADC64, 16 ns/sample).  Two trees,
// `normal` and `pedestal` (0.1 Hz random trigger), names as in produce_fadc.cc.
// Algorithm, settings and output: README.md.
//
//   produce_iadc -c prod_iadc.yaml -o out/prod_IADC_<run>.<sub>.root
//                /data/amore2/data/RAW/<run>/IADC_<run>.root.<sub>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <getopt.h>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TNamed.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"

#include "EventInfo.hh"
#include "FChannel.hh"
#include "FChannelData.hh"
#include "common.hh"
#include "yaml-cpp/yaml.h"

// FChannel waveforms, 12-bit, NDP = 8 RL - 2, read per channel.
// Zero suppression: physics events carry only the channels that fired.
constexpr int kMaxNdp = 4096; // buffer size: longer records are cut here
constexpr int kMaxCh = 520;

// Built-in defaults; prod_iadc.yaml overrides them.
constexpr int kBaseLoDef = 0, kBaseHiDef = 14; // baseline region, normal triggers
constexpr int kIntLoDef = 16, kIntHiDef = 48;
constexpr double kCF = 0.5; // constant fraction for the leading-edge time
constexpr double kSampleNs = 16.0; // IADC sample period
constexpr int kAdcCeiling = 4095;  // 12-bit ADC: samples here are clipped

// Charge level: header pedestal + qbase, qbase from the dark-pulse fit over
// window_normal (normal.baseline: fit; 0 if it fails) or 0 (header).
struct Config {
  int normLo = kBaseLoDef, normHi = kBaseHiDef; // normal: region the level comes from
  int pedLo = 0, pedHi = 0;                     // pedestal: same, 0 for the whole record
  int intLo = kIntLoDef, intHi = kIntHiDef;
  double cf = kCF;

  // dark-pulse fit on the pedestal triggers, see DarkFit below
  bool fit = true;
  // baseline under the charge: the dark-pulse fit over window_normal ("fit"),
  // or the ADC header pedestal ("header")
  bool normFit = true;
  double sigE = 2.0, seedH = 12.0, addRes = 8.0, addDchi2 = 25.0;
  int maxPulses = 12;
  double sigLo = 1.4, sigHi = 2.4, muMin = -8.0, keepH = 5.7;

  TString path = "(built-in defaults)";
};

bool LoadConfig(const char * path, Config & c)
{
  if (!path || !path[0]) return true;
  YAML::Node n;
  if (!LoadYaml(path, n)) return false;
  const auto pair = ReadPair;
  if (const YAML::Node b = n["baseline"]) {
    pair(b["window_normal"], c.normLo, c.normHi);
    pair(b["window_pedestal"], c.pedLo, c.pedHi);
  }
  if (const YAML::Node m = n["normal"]) {
    pair(m["integration"], c.intLo, c.intHi);
    if (m["cf"]) c.cf = m["cf"].as<double>();
    if (m["baseline"]) {
      const std::string v = m["baseline"].as<std::string>();
      if (v != "fit" && v != "header") {
        printf("ERROR: normal.baseline must be fit or header, not '%s'\n", v.c_str());
        return false;
      }
      c.normFit = v == "fit";
    }
  }
  if (const YAML::Node f = n["pedestal"] ? n["pedestal"]["fit"] : YAML::Node()) {
    if (f["enabled"]) c.fit = f["enabled"].as<bool>();
    if (f["sigma_e"]) c.sigE = f["sigma_e"].as<double>();
    if (f["seed_height"]) c.seedH = f["seed_height"].as<double>();
    if (f["add_residual"]) c.addRes = f["add_residual"].as<double>();
    if (f["add_dchi2"]) c.addDchi2 = f["add_dchi2"].as<double>();
    if (f["max_pulses"]) c.maxPulses = f["max_pulses"].as<int>();
    if (f["keep_height"]) c.keepH = f["keep_height"].as<double>();
    if (f["mu_min"]) c.muMin = f["mu_min"].as<double>();
    if (const YAML::Node v = f["sigma_range"]) {
      if (v.IsSequence() && v.size() == 2) {
        c.sigLo = v[0].as<double>();
        c.sigHi = v[1].as<double>();
      }
    }
  }
  c.path = path;
  return true;
}

// What both trees measure on a channel; no even/odd split (the IADC has none).
struct ChanResult {
  float base, rms, qpre;          // baseline region
  float hmax, tlead, tcf, charge; // integration window
  short hmaxx;
  short nsat; // samples of the record at the ADC ceiling
};

ChanResult Measure(const unsigned short * w, int nd, double ped, int lo, int hi, const Config & c)
{
  ChanResult r{};
  for (int i = 0; i < nd; i++)
    if (w[i] >= kAdcCeiling) r.nsat++;
  double s = 0, s2 = 0;
  for (int i = lo; i < hi; i++) {
    s += w[i];
    s2 += double(w[i]) * w[i];
  }
  const int np = hi - lo;
  const double mb = s / np;
  r.base = mb;
  r.rms = std::sqrt(std::max(0.0, s2 / np - mb * mb));
  r.qpre = 0;
  for (int i = c.normLo; i < c.normHi; i++)
    r.qpre += w[i] - ped;

  double q = 0, hm = -1e9;
  int hx = c.intLo;
  for (int i = c.intLo; i < c.intHi; i++) {
    const double d = w[i] - ped;
    q += d;
    if (d > hm) {
      hm = d;
      hx = i;
    }
  }
  r.charge = q;
  r.hmax = hm;
  r.hmaxx = hx;

  // first crossing of cf x hmax on the way up to the peak, interpolated
  r.tlead = -1;
  if (hm > 0) {
    const double lvl = c.cf * hm;
    for (int i = c.intLo + 1; i <= hx; i++) {
      const double a = w[i - 1] - ped, b = w[i] - ped;
      if (a < lvl && b >= lvl) {
        r.tlead = (i - 1) + (lvl - a) / (b - a);
        break;
      }
    }
  }

  // tcf: the cf crossing walking down from the peak (an earlier, smaller pulse
  // cannot take it), allowed before the window, in ns; not when saturated.
  r.tcf = -1;
  if (hm > 0 && r.nsat == 0) {
    const double lvl = c.cf * hm;
    int k = hx;
    while (k > 0 && w[k - 1] - ped > lvl)
      k--;
    if (k > 0) {
      const double a = w[k - 1] - ped, b = w[k] - ped;
      r.tcf = kSampleNs * ((k - 1) + (lvl - a) / (b - a));
    }
  }
  return r;
}

// ── Dark-pulse fit ─────────────────────────────────────────────────────────
// Flat level + n Moyal pulses of a common width,
//   b + sum_i H_i exp(-(z_i + exp(-z_i))/2 + 1/2),    z_i = (t - mu_i) / sigma
// area 4.1327 H sigma.  Seeds from find_peaks, then one pulse at a time where
// the residual is largest if it buys addDchi2; Levenberg-Marquardt with analytic
// derivatives, parameters clamped to their bounds.
constexpr double kMoyalArea = 4.1327; // sqrt(2 pi) e^(1/2)

struct DarkPulse {
  double h, mu, area;
};
struct DarkFitResult {
  bool ok = false;
  double b = 0, sigma = 0, chi2ndf = 0;
  std::vector<DarkPulse> pulses;
};

class DarkFit {
public:
  explicit DarkFit(const Config & c)
    : fC(c)
  {
  }

  // Fit samples [lo, hi) of the record y; positions come back in record samples.
  DarkFitResult Fit(const double * y, int lo, int hi) const
  {
    DarkFitResult out;
    const double * w = y + lo;
    const int n = hi - lo, t0 = lo;
    if (n < 4) return out;
    const double b0 = ClippedLevel(w, n);
    std::vector<double> r(n);
    for (int i = 0; i < n; i++)
      r[i] = w[i] - b0;

    std::vector<std::pair<double, double>> seeds; // (height, position)
    for (int i : FindPeaks(r))
      seeds.emplace_back(r[i], t0 + i);
    if (r[0] > fC.seedH && r[0] >= r[1]) seeds.emplace_back(1.5 * r[0], t0 - 1.0);
    if (r[n - 1] > fC.seedH && r[n - 1] >= r[n - 2]) seeds.emplace_back(r[n - 1], double(t0 + n));

    std::vector<double> p;
    double chi2 = Run(w, n, t0, b0, seeds, p);
    if (!std::isfinite(chi2)) return out;
    while ((int)seeds.size() < fC.maxPulses) {
      int imax = 0;
      double rmax = -1e9;
      for (int i = 0; i < n; i++) {
        const double d = w[i] - Model(t0 + i, p);
        if (d > rmax) {
          rmax = d;
          imax = i;
        }
      }
      if (rmax < fC.addRes) break;
      std::vector<std::pair<double, double>> s2;
      for (size_t j = 0; j < seeds.size(); j++)
        s2.emplace_back(p[2 + 2 * j], p[3 + 2 * j]);
      s2.emplace_back(rmax, double(t0 + imax));
      std::vector<double> p2;
      const double c2 = Run(w, n, t0, b0, s2, p2);
      if (!std::isfinite(c2) || chi2 - c2 < fC.addDchi2) break;
      seeds = s2;
      p = p2;
      chi2 = c2;
    }

    out.ok = true;
    out.b = p[0];
    out.sigma = p[1];
    out.chi2ndf = chi2 / std::max(1, n - (int)p.size());
    for (size_t j = 0; j < seeds.size(); j++) {
      const double h = p[2 + 2 * j], mu = p[3 + 2 * j];
      if (h > fC.keepH) out.pulses.push_back({h, mu, kMoyalArea * h * p[1]});
    }
    std::sort(out.pulses.begin(), out.pulses.end(),
              [](const DarkPulse & x, const DarkPulse & y) { return x.mu < y.mu; });
    return out;
  }

private:
  const Config & fC;

  // level from below: dark pulses only push samples up
  double ClippedLevel(const double * w, int n) const
  {
    std::vector<double> s(w, w + n);
    std::sort(s.begin(), s.end());
    const double pos = 0.2 * (n - 1);
    const int k = (int)pos;
    double m = s[k] + (pos - k) * (s[std::min(k + 1, n - 1)] - s[k]);
    for (int it = 0; it < 6; it++) {
      double sum = 0;
      int cnt = 0;
      for (int i = 0; i < n; i++)
        if (std::fabs(w[i] - m) < 2.5 * fC.sigE) {
          sum += w[i];
          cnt++;
        }
      if (cnt) m = sum / cnt;
    }
    return m;
  }

  std::vector<int> FindPeaks(const std::vector<double> & r) const
  {
    const int n = r.size();
    std::vector<int> cand;
    for (int i = 1; i < n - 1; i++)
      if (r[i] > r[i - 1] && r[i] >= r[i + 1] && r[i] >= fC.seedH) cand.push_back(i);
    // distance 3: the higher peak wins
    std::vector<int> order(cand);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return r[a] > r[b]; });
    std::vector<int> kept;
    for (int i : order) {
      bool close = false;
      for (int k : kept)
        if (std::abs(i - k) < 3) close = true;
      if (!close) kept.push_back(i);
    }
    std::vector<int> out;
    for (int i : kept) {
      double lmin = r[i], rmin = r[i];
      for (int k = i - 1; k >= 0 && r[k] <= r[i]; k--)
        lmin = std::min(lmin, r[k]);
      for (int k = i + 1; k < n && r[k] <= r[i]; k++)
        rmin = std::min(rmin, r[k]);
      if (r[i] - std::max(lmin, rmin) >= 0.8 * fC.seedH) out.push_back(i);
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  static double Model(double t, const std::vector<double> & p)
  {
    double y = p[0];
    for (size_t j = 2; j + 1 < p.size(); j += 2) {
      const double z = std::max((t - p[j + 1]) / p[1], -30.0);
      y += p[j] * std::exp(-0.5 * (z + std::exp(-z)) + 0.5);
    }
    return y;
  }

  double Chi2(const double * w, int n, int t0, const std::vector<double> & p) const
  {
    double c = 0;
    for (int i = 0; i < n; i++) {
      const double d = (w[i] - Model(t0 + i, p)) / fC.sigE;
      c += d * d;
    }
    return c;
  }

  // Levenberg-Marquardt from these seeds; returns chi2 and the parameters
  double Run(const double * w, int n, int t0, double b0,
             const std::vector<std::pair<double, double>> & seeds, std::vector<double> & p) const
  {
    const int np = 2 + 2 * seeds.size();
    std::vector<double> lo(np), hi(np);
    lo[0] = b0 - 30;
    hi[0] = b0 + 30;
    lo[1] = fC.sigLo;
    hi[1] = fC.sigHi;
    p.assign(np, 0);
    p[0] = b0;
    p[1] = 1.8;
    for (size_t j = 0; j < seeds.size(); j++) {
      lo[2 + 2 * j] = 0;
      hi[2 + 2 * j] = 5000;
      lo[3 + 2 * j] = t0 + fC.muMin;
      hi[3 + 2 * j] = t0 + n + 10;
      p[2 + 2 * j] = seeds[j].first;
      p[3 + 2 * j] = seeds[j].second;
    }
    auto clamp = [&](std::vector<double> & q) {
      for (int k = 0; k < np; k++)
        q[k] = std::min(std::max(q[k], lo[k] + 1e-6), hi[k] - 1e-6);
    };
    clamp(p);

    double chi2 = Chi2(w, n, t0, p), lambda = 1e-3;
    std::vector<double> J(n * np), A(np * np), g(np), dp(np), q(np);
    for (int it = 0; it < 300; it++) {
      // residuals and analytic Jacobian
      std::fill(A.begin(), A.end(), 0.0);
      std::fill(g.begin(), g.end(), 0.0);
      for (int i = 0; i < n; i++) {
        double * Ji = &J[i * np];
        Ji[0] = 1;
        Ji[1] = 0;
        double m = p[0];
        for (int j = 2; j < np; j += 2) {
          const double zr = (t0 + i - p[j + 1]) / p[1];
          const double z = std::max(zr, -30.0);
          const double ez = std::exp(-z);
          const double E = std::exp(-0.5 * (z + ez) + 0.5);
          m += p[j] * E;
          const double dEdz = zr > -30.0 ? -0.5 * (1 - ez) * E : 0.0;
          Ji[j] = E;
          Ji[j + 1] = -p[j] * dEdz / p[1];
          Ji[1] += -p[j] * dEdz * z / p[1];
        }
        const double ri = w[i] - m;
        for (int a = 0; a < np; a++) {
          g[a] += Ji[a] * ri;
          for (int b = 0; b <= a; b++)
            A[a * np + b] += Ji[a] * Ji[b];
        }
      }
      for (int a = 0; a < np; a++)
        for (int b = 0; b < a; b++)
          A[b * np + a] = A[a * np + b];

      bool accepted = false;
      double c2n = chi2;
      while (lambda < 1e10) {
        if (!Solve(A, g, lambda, np, dp)) {
          lambda *= 10;
          continue;
        }
        for (int k = 0; k < np; k++)
          q[k] = p[k] + dp[k];
        clamp(q);
        c2n = Chi2(w, n, t0, q);
        if (c2n < chi2) {
          accepted = true;
          lambda = std::max(lambda / 10, 1e-9);
          break;
        }
        lambda *= 10;
      }
      if (!accepted) break;
      const double gain = chi2 - c2n;
      p = q;
      chi2 = c2n;
      if (gain < 1e-6 * chi2 + 1e-9) break;
    }
    return chi2;
  }

  // (A + lambda diag A) x = g, Gaussian elimination with partial pivoting
  static bool Solve(const std::vector<double> & A, const std::vector<double> & g, double lambda,
                    int np, std::vector<double> & x)
  {
    std::vector<double> M(np * (np + 1));
    for (int a = 0; a < np; a++) {
      for (int b = 0; b < np; b++)
        M[a * (np + 1) + b] = A[a * np + b];
      M[a * (np + 1) + a] += lambda * (A[a * np + a] + 1e-12);
      M[a * (np + 1) + np] = g[a];
    }
    for (int c = 0; c < np; c++) {
      int piv = c;
      for (int r2 = c + 1; r2 < np; r2++)
        if (std::fabs(M[r2 * (np + 1) + c]) > std::fabs(M[piv * (np + 1) + c])) piv = r2;
      if (std::fabs(M[piv * (np + 1) + c]) < 1e-300) return false;
      if (piv != c)
        for (int k = 0; k <= np; k++)
          std::swap(M[c * (np + 1) + k], M[piv * (np + 1) + k]);
      for (int r2 = c + 1; r2 < np; r2++) {
        const double f = M[r2 * (np + 1) + c] / M[c * (np + 1) + c];
        for (int k = c; k <= np; k++)
          M[r2 * (np + 1) + k] -= f * M[c * (np + 1) + k];
      }
    }
    x.assign(np, 0);
    for (int c = np - 1; c >= 0; c--) {
      double sum = M[c * (np + 1) + np];
      for (int k = c + 1; k < np; k++)
        sum -= M[c * (np + 1) + k] * x[k];
      x[c] = sum / M[c * (np + 1) + c];
    }
    return true;
  }
};

int produce(const char * in, const char * out, const char * cfgPath)
{
  TString stream;
  int run = -1, sub = -1;
  if (!ParseRawName(in, stream, run, sub) || stream != "IADC") {
    printf("ERROR: '%s' is not named IADC_<run>.root.<subrun>\n", in);
    return 1;
  }
  Config cfg;
  if (!LoadConfig(cfgPath, cfg)) return 1;
  // pedHi 0 = the whole record, whatever its length on the channel
  const int pedLo = cfg.pedLo, pedHi = cfg.pedHi;
  if (cfg.normHi <= cfg.normLo || cfg.intHi <= cfg.intLo || cfg.intHi > kMaxNdp ||
      (pedHi > 0 && pedHi <= pedLo) || pedHi > kMaxNdp) {
    printf("ERROR: bad windows: baseline normal [%d,%d) pedestal [%d,%d)  integration [%d,%d)\n",
           cfg.normLo, cfg.normHi, pedLo, pedHi, cfg.intLo, cfg.intHi);
    return 1;
  }

  TFile * fin = TFile::Open(in);
  if (!fin || fin->IsZombie()) {
    printf("ERROR: cannot open %s\n", in);
    return 1;
  }
  auto * tin = static_cast<TTree *>(fin->Get("AbsEvent"));
  if (!tin) {
    printf("ERROR: no AbsEvent tree in %s\n", in);
    return 1;
  }
  EventInfo * info = nullptr;
  FChannelData * data = nullptr;
  tin->SetBranchAddress("EventInfo", &info);
  tin->SetBranchAddress("FChannelData", &data);

  printf("run %d  subrun %d\n", run, sub);
  printf("config : %s\n", cfg.path.Data());
  printf("input  : %s  (%lld events)\n", in, tin->GetEntries());
  printf("output : %s\n", out);
  printf("normal   : baseline [%d,%d)  integration [%d,%d)  cf %.2f\n", cfg.normLo, cfg.normHi,
         cfg.intLo, cfg.intHi, cfg.cf);
  printf("normal   : charge baseline %s\n",
         cfg.normFit ? "dark-pulse fit over window_normal" : "header pedestal");
  if (pedHi > 0) printf("pedestal : baseline [%d,%d)\n", pedLo, pedHi);
  else printf("pedestal : baseline [%d, end of record)\n", pedLo);
  if (cfg.fit)
    printf("pedestal : dark-pulse fit, sigma_e %.1f, seed %.1f, add %.1f / dchi2 %.0f, "
           "<= %d pulses, width [%.1f,%.1f], keep h > %.1f\n",
           cfg.sigE, cfg.seedH, cfg.addRes, cfg.addDchi2, cfg.maxPulses, cfg.sigLo, cfg.sigHi,
           cfg.keepH);

  TFile * fout = TFile::Open(out, "recreate", "", 1);
  if (!fout || fout->IsZombie()) {
    printf("ERROR: cannot write %s\n", out);
    return 1;
  }

  // Shared quantities have the same name and type in both trees.
  Int_t runNo = run, subrun = sub, nch = 0, nhit = 0, nmod = 0;
  UInt_t evtNum = 0, trgNum = 0;
  ULong64_t trgTime = 0;
  Short_t pid[kMaxCh], hmaxx[kMaxCh], nsat[kMaxCh];
  UShort_t tbit[kMaxCh], ped[kMaxCh];
  Float_t base[kMaxCh], rms[kMaxCh], qpre[kMaxCh], hmax[kMaxCh], tlead[kMaxCh], tcf[kMaxCh];
  Float_t charge[kMaxCh];
  Float_t qbase[kMaxCh]; // ADC over ped: the level charge, hmax, tlead and qpre are taken against

  // the dark-pulse fit, per channel and per pulse: over window_normal on the
  // normal tree, over the whole record on the pedestal tree
  constexpr int kMaxPul = kMaxCh * 12;
  Short_t npul[kMaxCh];
  Float_t fit_base[kMaxCh], fit_sigma[kMaxCh], fit_chi2[kMaxCh];
  Int_t npulse = 0;
  Short_t pul_ch[kMaxPul];
  Float_t pul_t[kMaxPul], pul_h[kMaxPul], pul_area[kMaxPul];
  const DarkFit fitter(cfg);
  if (cfg.maxPulses > 12) {
    printf("ERROR: max_pulses %d above the 12 the arrays are sized for\n", cfg.maxPulses);
    return 1;
  }

  auto book = [&](TTree * t, bool isPed) {
    if (!isPed) t->Branch("run", &runNo, "run/I");
    t->Branch("nch", &nch, "nch/I");
    t->Branch("subrun", &subrun, "subrun/I");
    t->Branch("evtNum", &evtNum, "evtNum/i");
    t->Branch("trgNum", &trgNum, "trgNum/i");
    t->Branch("trgTime", &trgTime, "trgTime/l"); // ns since the run started
    if (!isPed) {
      t->Branch("nhit", &nhit, "nhit/I"); // channels with the trigger bit set
      t->Branch("nmod", &nmod, "nmod/I"); // modules among the channels, (pid-1)/4
    }
    t->Branch("pid", pid, "pid[nch]/S");          // channel id, as the DAQ writes it
    t->Branch("tbit", tbit, "tbit[nch]/s");       // DAQ per-channel discriminator
    t->Branch("ped", ped, "ped[nch]/s");          // ADC, pedestal from the ADC header
    t->Branch("base", base, "base[nch]/F");       // ADC, mean of the baseline region
    t->Branch("rms", rms, "rms[nch]/F");          // ADC, its rms
    t->Branch("qpre", qpre, "qpre[nch]/F");       // ADC x sample, window_normal over ped
    t->Branch("hmax", hmax, "hmax[nch]/F");       // ADC over ped, in the integration window
    t->Branch("hmaxx", hmaxx, "hmaxx[nch]/S");    // where it is
    t->Branch("tlead", tlead, "tlead[nch]/F");    // sample, cf crossing; -1 if none
    t->Branch("tcf", tcf, "tcf[nch]/F");          // ns, cf crossing walking down from the peak; -1 if none or saturated
    t->Branch("nsat", nsat, "nsat[nch]/S");       // samples at the ADC ceiling (4095); > 0 = saturated
    t->Branch("charge", charge, "charge[nch]/F"); // ADC x sample, integration window over qbase
    t->Branch("qbase", qbase, "qbase[nch]/F");    // ADC over ped, the level of the line above
    if (isPed ? cfg.fit : cfg.normFit) {
      t->Branch("fit_base", fit_base, "fit_base[nch]/F");    // ADC over ped, fitted flat level
      t->Branch("fit_sigma", fit_sigma, "fit_sigma[nch]/F"); // sample, common Moyal width
      t->Branch("fit_chi2", fit_chi2, "fit_chi2[nch]/F");    // chi2/ndf, white noise sigma_e
      t->Branch("npul", npul, "npul[nch]/S");                // pulses kept, -1 if the fit failed
      t->Branch("npulse", &npulse, "npulse/I");              // pulses in the event
      t->Branch("pul_ch", pul_ch, "pul_ch[npulse]/S");       // index into the [nch] arrays
      t->Branch("pul_t", pul_t, "pul_t[npulse]/F");          // sample, Moyal peak position
      t->Branch("pul_h", pul_h, "pul_h[npulse]/F");          // ADC, peak height
      t->Branch("pul_area", pul_area, "pul_area[npulse]/F"); // ADC x sample, 4.1327 h sigma
    }
  };
  auto * tnorm = new TTree("normal", "PSMD normal-trigger production");
  auto * tped = new TTree("pedestal", "PSMD pedestal-trigger production");
  book(tnorm, false);
  book(tped, true);

  TStopwatch watch;
  Long64_t nshort = 0;
  const Long64_t n = tin->GetEntries();
  for (Long64_t i = 0; i < n; i++) {
    tin->GetEntry(i);
    const bool isPed = info->GetTriggerType() == kPedTrigType;
    evtNum = info->GetEventNumber();
    trgNum = info->GetTriggerNumber();
    trgTime = info->GetTriggerTime();
    nch = 0;
    nhit = 0;
    npulse = 0;
    std::set<int> mods;
    for (int j = 0; j < data->GetN() && nch < kMaxCh; j++) {
      const FChannel * ch = data->Get(j);
      const int nd = std::min(ch->GetNdp(), kMaxNdp);
      const int pHi = pedHi > 0 ? pedHi : nd;
      if (nd < cfg.intHi || nd < cfg.normHi || nd < pHi) {
        nshort++;
        continue;
      }
      const unsigned short * wv = ch->GetWaveform();
      const double hped = ch->GetPedestal();
      double y[kMaxNdp];
      for (int k = 0; k < nd; k++)
        y[k] = wv[k] - hped;

      // one fit gives the level and the dark pulses
      double lvl = 0;
      DarkFitResult fn;
      if (isPed ? (cfg.fit || cfg.normFit) : cfg.normFit)
        fn = isPed ? fitter.Fit(y, 0, nd) : fitter.Fit(y, cfg.normLo, cfg.normHi);
      if (cfg.normFit && fn.ok) lvl = fn.b; // a failed fit falls back to the header pedestal
      const ChanResult r = isPed ? Measure(wv, nd, hped + lvl, pedLo, pHi, cfg)
                                 : Measure(wv, nd, hped + lvl, cfg.normLo, cfg.normHi, cfg);
      qbase[nch] = lvl;
      pid[nch] = ch->GetID();
      tbit[nch] = ch->GetBit();
      ped[nch] = ch->GetPedestal();
      base[nch] = r.base;
      rms[nch] = r.rms;
      qpre[nch] = r.qpre;
      hmax[nch] = r.hmax;
      hmaxx[nch] = r.hmaxx;
      tlead[nch] = r.tlead;
      tcf[nch] = r.tcf;
      nsat[nch] = r.nsat;
      charge[nch] = r.charge;
      if (isPed ? cfg.fit : cfg.normFit) {
        const DarkFitResult & f = fn;
        fit_base[nch] = f.ok ? f.b : 0;
        fit_sigma[nch] = f.ok ? f.sigma : 0;
        fit_chi2[nch] = f.ok ? f.chi2ndf : -1;
        npul[nch] = f.ok ? f.pulses.size() : -1;
        for (const DarkPulse & dp : f.pulses) {
          if (npulse >= kMaxPul) break;
          pul_ch[npulse] = nch;
          pul_t[npulse] = dp.mu;
          pul_h[npulse] = dp.h;
          pul_area[npulse] = dp.area;
          npulse++;
        }
      }
      mods.insert((ch->GetID() - 1) / 4);
      if (ch->GetBit() != 0) nhit++;
      nch++;
    }
    nmod = mods.size();
    (isPed ? tped : tnorm)->Fill();
  }

  fout->cd();
  tnorm->Write();
  tped->Write();
  TNamed("run", Form("%d", run)).Write();
  TNamed("subrun", Form("%d", sub)).Write();
  TNamed("input", in).Write();
  TNamed("config", cfg.path.Data()).Write();
  TNamed("window_normal", Form("[%d,%d)", cfg.normLo, cfg.normHi)).Write();
  TNamed("window_pedestal", pedHi > 0 ? Form("[%d,%d)", pedLo, pedHi) : Form("[%d,end)", pedLo))
      .Write();
  TNamed("integration", Form("[%d,%d)", cfg.intLo, cfg.intHi)).Write();
  TNamed("cf", Form("%.3f", cfg.cf)).Write();
  TNamed("baseline", cfg.normFit ? "header pedestal + dark-pulse fit over window_normal (qbase)"
                                 : "header pedestal")
      .Write();
  TNamed("sample_ns", Form("%g", kSampleNs)).Write();
  TNamed("adc_ceiling", Form("%d", kAdcCeiling)).Write();
  TNamed("pedestal_fit", cfg.fit
                             ? Form("moyal+flat: sigma_e %.2f seed %.2f add %.2f dchi2 %.1f max %d "
                                    "width [%.2f,%.2f] mu_min %.1f keep_h %.2f",
                                    cfg.sigE, cfg.seedH, cfg.addRes, cfg.addDchi2, cfg.maxPulses,
                                    cfg.sigLo, cfg.sigHi, cfg.muMin, cfg.keepH)
                             : "off")
      .Write();
  if (nshort > 0) printf("WARNING: %lld waveforms were too short and were skipped\n", nshort);
  printf("done: %lld events read -> %lld normal + %lld pedestal, in %.1f s\n", n,
         tnorm->GetEntries(), tped->GetEntries(), watch.RealTime());
  fout->Close();
  fin->Close();
  return 0;
}

// ── Command line ───────────────────────────────────────────────────────────
static void Usage(const char * prog)
{
  printf("usage: %s -o OUTPUT [-c CONFIG] INPUT\n"
         "  INPUT               raw IADC_<run>.root.<subrun> file\n"
         "  -o, --out FILE      output file (required)\n"
         "  -c, --config FILE   YAML settings (default: built-in)\n",
         prog);
}

int main(int argc, char ** argv)
{
  const char *cfg = "", *out = nullptr;

  const option opts[] = {{"out", required_argument, nullptr, 'o'},
                         {"config", required_argument, nullptr, 'c'},
                         {"help", no_argument, nullptr, 'h'},
                         {nullptr, 0, nullptr, 0}};
  int c;
  while ((c = getopt_long(argc, argv, "o:c:h", opts, nullptr)) != -1) {
    switch (c) {
      case 'o': out = optarg; break;
      case 'c': cfg = optarg; break;
      case 'h': Usage(argv[0]); return 0;
      default: Usage(argv[0]); return 2;
    }
  }
  if (!out || argc - optind != 1) {
    Usage(argv[0]);
    return 2;
  }
  return produce(argv[optind], out, cfg);
}
