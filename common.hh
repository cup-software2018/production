// What produce_fadc and produce_iadc share: which trigger is the pedestal one,
// reading the DAQ's file names, the YAML settings and the command line.
#pragma once

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <exception>

#include "TString.h"
#include "TSystem.h"

#include "yaml-cpp/yaml.h"

// trgType 1 is the 0.1 Hz pedestal trigger on every board; everything else is
// a normal trigger.
constexpr unsigned short kPedTrigType = 1;

// The DAQ names its files <stream>_<run:06d>.root.<subrun:05d>, e.g.
// FADC1_000123.root.00001 or IADC_000123.root.00001; the run and subrun the
// output carries are read from there.
inline bool ParseRawName(const char * path, TString & stream, int & run, int & sub)
{
  const TString name = gSystem->BaseName(path);
  const Ssiz_t us = name.Index("_");
  if (us <= 0) return false;
  stream = name(0, us);
  int n = 0;
  return sscanf(name.Data() + us + 1, "%d.root.%d%n", &run, &sub, &n) == 2 &&
         n == name.Length() - us - 1;
}

// Parse the settings file; false, with the reason printed, if it is missing
// or broken.
inline bool LoadYaml(const char * path, YAML::Node & n)
{
  if (gSystem->AccessPathName(path)) {
    printf("ERROR: no config file '%s'\n", path);
    return false;
  }
  try {
    n = YAML::LoadFile(path);
  } catch (const std::exception & e) {
    printf("ERROR: cannot parse '%s': %s\n", path, e.what());
    return false;
  }
  return true;
}

// A two-element sequence [lo, hi]; anything else leaves lo and hi alone, so a
// partial config is a set of overrides rather than a trap.
inline void ReadPair(const YAML::Node & v, int & lo, int & hi)
{
  if (v && v.IsSequence() && v.size() == 2) {
    lo = v[0].as<int>();
    hi = v[1].as<int>();
  }
}

// A whole decimal number from the command line, or false.
inline bool ParseInt(const char * s, long long & v)
{
  if (!s || !s[0]) return false;
  char * end = nullptr;
  errno = 0;
  v = std::strtoll(s, &end, 10);
  return errno == 0 && *end == '\0';
}
