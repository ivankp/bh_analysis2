#ifndef BH_ANALYSIS_PARSE_BINS_FILE_NAME_HH
#define BH_ANALYSIS_PARSE_BINS_FILE_NAME_HH

#include <cstring>

namespace parse {

bool bins_file_name(const char* arg, const char*& file) {
  const size_t len = std::strlen(arg);

  if (len<6) return false;
  if (std::strcmp(arg+(len-5),".bins")) return false;

  if (file) {
    std::cout << "\n\033[31mWarning: "
                 "More than 1 argument interpreted as binning file name\n"
              << file << "\nreplaced by\n" << (arg+4)
              << "\n\033[0m" << std::endl;
  }
  
  file = arg+4;
  return true;
}

}

#endif
