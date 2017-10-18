#ifndef BH_ANALYSIS_PARSE_JET_PT_CUT_HH
#define BH_ANALYSIS_PARSE_JET_PT_CUT_HH

#include <cstdlib>
#include <cctype>

namespace parse {

bool jet_pt_cut(const char* arg, double& cut) {
  const char * const str = "jetpt";

  char c = '\0';
  int i = 0;
  for ( ; (c = arg[i]) && i<5; ++i)
    if (c != str[i]) return false;
  if (i!=5) return false;
  
  cut = std::atof(arg+i);
  return true;
}

}

#endif

