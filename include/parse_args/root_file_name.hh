#ifndef BH_ANALYSIS_PARSE_ROOT_FILE_NAME_HH
#define BH_ANALYSIS_PARSE_ROOT_FILE_NAME_HH

#include <cstring>

namespace parse {

bool root_file_name(const char* arg, std::vector<const char*>& files) {
  const size_t len = std::strlen(arg);

  if (len<6) return false;
  if (std::strcmp(arg+(len-5),".root")) return false;
  
  files.push_back(arg);
  return true;
}

}

#endif
