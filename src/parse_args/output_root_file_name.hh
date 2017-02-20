#ifndef BH_ANALYSIS_PARSE_OUTPUT_ROOT_FILE_NAME_HH
#define BH_ANALYSIS_PARSE_OUTPUT_ROOT_FILE_NAME_HH

#include <cstring>

namespace parse {

bool output_root_file_name(const char* arg, const char*& file) {
  const size_t len = std::strlen(arg);

  if (len<5) return false;
  for (int i=0; i<4; ++i)
    if (arg[i] != "out:"[i]) return false;
  
  if (std::strcmp(arg+(len-5),".root")) {
    if (std::strcmp(arg+4,"/dev/null")) { // allow output to /dev/null
      std::cerr << "\n\033[31mOutput file must have .root extension\033[0m"
                << std::endl;
      return false;
    }
  }

  if (file) {
    std::cout << "\n\033[31mWarning: "
                 "More than 1 argument interpreted as output file name\n"
              << file << "\nreplaced by\n" << (arg+4)
              << "\033[0m\n" << std::endl;
  }
  
  file = arg+4;
  return true;
}

}

#endif
