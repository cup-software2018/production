// What produce_fadc and produce_iadc share: which trigger is the pedestal one,
// reading the DAQ's file names, the YAML settings and the command line.
#pragma once

#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <sstream>
#include <string>

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
// or broken.  $NAME or ${NAME} in it is replaced by the environment variable
// NAME, which must be set (prod_wcmd_subrun.sh and prod_psmd_subrun.sh set
// the windows this way).
inline bool LoadYaml(const char * path, YAML::Node & n)
{
  if (gSystem->AccessPathName(path)) {
    printf("ERROR: no config file '%s'\n", path);
    return false;
  }
  std::ifstream f(path);
  std::stringstream ss;
  ss << f.rdbuf();
  std::string text = ss.str();
  const auto isName = [](char ch, bool first) {
    return ch == '_' || std::isalpha((unsigned char)ch) || (!first && std::isdigit((unsigned char)ch));
  };
  for (size_t at = text.find('$'); at != std::string::npos; at = text.find('$', at)) {
    size_t from, len, skip; // name at [from, from + len); the whole reference is skip long
    if (at + 1 < text.size() && text[at + 1] == '{') {
      const size_t end = text.find('}', at);
      if (end == std::string::npos) break;
      from = at + 2, len = end - from, skip = end + 1 - at;
    } else if (at + 1 < text.size() && isName(text[at + 1], true)) {
      size_t end = at + 2;
      while (end < text.size() && isName(text[end], false))
        end++;
      from = at + 1, len = end - from, skip = end - at;
    } else {
      at++; // a lone $
      continue;
    }
    const std::string name = text.substr(from, len);
    const char * v = std::getenv(name.c_str());
    if (!v) {
      printf("ERROR: '%s' uses $%s, which is not set\n", path, name.c_str());
      return false;
    }
    text.replace(at, skip, v);
    at += std::strlen(v);
  }
  try {
    n = YAML::Load(text);
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
