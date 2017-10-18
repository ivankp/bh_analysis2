#ifndef BH_ANALYSIS_PARSE_TREE_NAME_HH
#define BH_ANALYSIS_PARSE_TREE_NAME_HH

#include <cstring>

namespace parse {

bool tree_name(const char* arg, const char*& name) {
  for (int i=0; i<5; ++i)
    if (arg[i] != "tree:"[i]) return false;
  if (arg[5]=='\n') return false;
  
  name = arg+5;
  return true;
}

}

#endif
