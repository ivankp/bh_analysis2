#ifndef BH_ANALYSIS_PARSE_OUTPUT_DAT_FILE_NAME_HH
#define BH_ANALYSIS_PARSE_OUTPUT_DAT_FILE_NAME_HH

#include <cstring>

namespace parse {

bool output_dat_file_name(const char* arg, const char*& file) {
  const size_t len = std::strlen(arg);

  if (len<5) return false;
  
  if (std::strcmp(arg+(len-4),".dat")) return false;

  if (file) {
    std::cout << "\n\033[31mWarning: "
                 "More than 1 argument interpreted as output file name\n"
              << file << "\nreplaced by\n" << arg
              << "\033[0m\n" << std::endl;
  }
  
  file = arg;
  return true;
}

}

#endif

