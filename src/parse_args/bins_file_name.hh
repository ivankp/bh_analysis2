#ifndef BH_ANALYSIS_PARSE_BINS_FILE_NAME_HH
#define BH_ANALYSIS_PARSE_BINS_FILE_NAME_HH

#include <cstring>

namespace parse {

bool bins_file_name(const char* arg, const char*& file) {
  const size_t len = std::strlen(arg);

  if (len<6) return false;
  if (std::strcmp(arg+(len-5),".bins")) return false;

  file = arg;
  return true;
}

}

#endif
