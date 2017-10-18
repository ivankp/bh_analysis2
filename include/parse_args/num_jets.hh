#ifndef BH_ANALYSIS_PARSE_NUM_JETS_HH
#define BH_ANALYSIS_PARSE_NUM_JETS_HH

#include <cstdlib>
#include <cctype>

namespace parse {

bool num_jets(const char* arg, unsigned& njets) {
  char c = '\0';
  int i = 0;
  for ( ; (c = arg[i]); ++i)
    if (c=='\0' || !std::isdigit(c)) break;

  if (c!='j') return false;
  ++i;
  if (arg[i]!='\0') return false;

  const auto prev_njets = njets;
  njets = std::atoi(arg);

  if (prev_njets) {
    std::cout << "\n\033[31mWarning: "
                 "More than 1 argument interpreted as number of jets\n"
              << prev_njets << "j replaced by "
              << njets << "j\033[0m\n" << std::endl;
  }

  return true;
}

}

#endif

