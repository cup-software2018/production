// Read the PSMD gain-matching table (psmd_gain_<run>.txt, from 3_muon_cal.py) and apply it;
// the C++ / ROOT-macro counterpart of psmd_gain.py.
//
//     q = charge / (g[pid-1] * c[pid-1])
//
//   #include "psmd_gain.hh"
//   psmd_gain::Table t = psmd_gain::Load("psmd_gain_<run>.txt");   // by PID-1: g, c, use
//   double q = t.Calibrate(charge, pid);                         // NaN where use == 0
//   auto d = psmd_gain::Define(df, t);                           // RDataFrame: q[nch] next to
//   charge[nch], 0 where use == 0 ROOT::RDataFrame df("normal", psmd_gain::ProdFiles(prod, run));
//
// Header only; works in a ROOT macro (.L / #include) or compiled code.
#pragma once

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "TString.h"
#include "TSystem.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

namespace psmd_gain
{

constexpr int kNPid = 520;

struct Table {
  std::vector<double> g, c; // by PID-1
  std::vector<bool> use;

  // g * c, NaN where the channel is not usable (or pid out of range)
  double Scale(int pid) const
  {
    if (pid < 1 || pid > static_cast<int>(use.size()) || !use[pid - 1])
      return std::numeric_limits<double>::quiet_NaN();
    return g[pid - 1] * c[pid - 1];
  }
  double Calibrate(double charge, int pid) const { return charge / Scale(pid); }

  // one event: q[nch], 0 where use == 0
  template <typename P, typename Q>
  ROOT::RVec<float> Calibrate(const ROOT::RVec<P> & pid, const ROOT::RVec<Q> & charge) const
  {
    ROOT::RVec<float> o(pid.size());
    for (size_t i = 0; i < pid.size(); i++) {
      double s = Scale(pid[i]);
      o[i] = std::isfinite(s) ? charge[i] / s : 0.f;
    }
    return o;
  }
};

// Glob of the IADC production files of a run: <prod>/<run:06d>/ (as PROD), or <prod>/ directly.
inline std::string ProdFiles(const std::string & prod, int run)
{
  TString sub = TString::Format("%s/%06d", prod.c_str(), run);
  Long_t id, flags, modtime;
  Long64_t size;
  bool isdir = gSystem->GetPathInfo(sub, &id, &size, &flags, &modtime) == 0 && (flags & 2);
  return TString::Format("%s/prod_IADC_%06d.*.root", isdir ? sub.Data() : prod.c_str(), run).Data();
}

inline Table Load(const std::string & path)
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("psmd_gain::Load: cannot open " + path);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Table t;
  t.g.assign(kNPid, nan);
  t.c.assign(kNPid, nan);
  t.use.assign(kNPid, false);
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::string spid, sg, sc, suse; // strings: "nan" does not parse as a double with >>
    if (!(ss >> spid >> sg >> sc >> suse)) continue;
    int pid = std::stoi(spid);
    if (pid < 1 || pid > kNPid) continue;
    t.g[pid - 1] = std::stod(sg);
    t.c[pid - 1] = std::stod(sc);
    t.use[pid - 1] = std::stoi(suse) != 0;
  }
  return t;
}

// q[nch] (or col) from the pid and charge columns
template <typename D>
auto Define(D df, const Table & t, const std::string & col = "q")
{
  return df.Define(col,
                   [t](const ROOT::RVec<short> & pid, const ROOT::RVec<float> & charge) {
                     return t.Calibrate(pid, charge);
                   },
                   {"pid", "charge"});
}

} // namespace psmd_gain
