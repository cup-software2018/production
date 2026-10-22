// Read the WCMD gain table (wcmd_gain_<run>.txt, from 2_spe_fit.py) and apply it;
// the C++ / ROOT-macro counterpart of wcmd_gain.py.
//
//     q = charge / g[pid-1]                 in photoelectrons
//
//   #include "wcmd_gain.hh"
//   wcmd_gain::Table t = wcmd_gain::Load("wcmd_gain_<run>.txt");   // by PID-1: g, use
//   double q = t.Calibrate(charge, pid);                         // NaN where use == 0
//   auto d = wcmd_gain::Define(df, t);                           // RDataFrame: q[nch] next to
//   charge[nch], 0 where use == 0 ROOT::RDataFrame df("normal", wcmd_gain::ProdFiles(prod, run));
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

namespace wcmd_gain
{

struct Table {
  std::vector<double> g; // by PID-1, up to the largest PID in the table
  std::vector<bool> use;

  // g, NaN where the channel is not usable (or pid out of range)
  double Scale(int pid) const
  {
    if (pid < 1 || pid > static_cast<int>(use.size()) || !use[pid - 1])
      return std::numeric_limits<double>::quiet_NaN();
    return g[pid - 1];
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

// Glob of the FADC production files of a run: <prod>/<run:06d>/ (as PROD), or <prod>/ directly.
inline std::string ProdFiles(const std::string & prod, int run, int board = 1)
{
  TString sub = TString::Format("%s/%06d", prod.c_str(), run);
  Long_t id, flags, modtime;
  Long64_t size;
  bool isdir = gSystem->GetPathInfo(sub, &id, &size, &flags, &modtime) == 0 && (flags & 2);
  return TString::Format("%s/prod_FADC%d_%06d.*.root", isdir ? sub.Data() : prod.c_str(), board,
                         run)
      .Data();
}

// columns: pid g g_err s_spe s_p rate bias chi2 ndf use
inline Table Load(const std::string & path)
{
  std::ifstream in(path);
  if (!in) throw std::runtime_error("wcmd_gain::Load: cannot open " + path);
  const double nan = std::numeric_limits<double>::quiet_NaN();
  Table t;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ss(line);
    std::vector<std::string> col; // strings: "nan" does not parse as a double with >>
    for (std::string w; ss >> w;)
      col.push_back(w);
    if (col.size() < 3) continue;
    int pid = std::stoi(col.front());
    if (pid < 1) continue;
    if (pid > static_cast<int>(t.use.size())) {
      t.g.resize(pid, nan);
      t.use.resize(pid, false);
    }
    t.g[pid - 1] = std::stod(col[1]);
    t.use[pid - 1] = std::stoi(col.back()) != 0;
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

} // namespace wcmd_gain
