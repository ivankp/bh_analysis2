#ifndef IVANP_PARSE_JETDEF_HH
#define IVANP_PARSE_JETDEF_HH

#include <cstdlib>
#include <cctype>

namespace parse {

bool jetdef(const char* arg,
#ifdef __FASTJET_JETDEFINITION_HH__
  fastjet
#elif defined(__FJCORE_JETDEFINITION_HH__)
  fjcore
#else
#error "JetDefinition requires fastjet or fjcore headers included"
#endif
  ::JetDefinition& def
) {
  char c;
  int i = 0;
  for ( ; (c = arg[i]); ++i)
    if (c=='\0' || !std::isalpha(c)) break;

  auto cmp = [arg,i](const char* str, int len) {
    if (i!=len) return false;
    for (int j=0; j<len; ++j)
      if (std::tolower(arg[j]) != str[j]) return false;
    return true;
  };

  using namespace
#ifdef __FASTJET_JETDEFINITION_HH__
  fastjet;
#else
  fjcore;
#endif

  JetAlgorithm alg;

  if (!std::isdigit(c)) return false;
  if (cmp("kt",2)) alg = kt_algorithm;
  else if (cmp("antikt",6) || cmp("akt",3)) alg = antikt_algorithm;
  else if (cmp("cambridge",9) || cmp("ca",2)) alg = cambridge_algorithm;
  else return false;

  const int r = i;
  for ( ; (c = arg[i]); ++i) {
    if (c=='\0') break;
    if (!std::isdigit(c)) return false;
  }

  const JetAlgorithm prev_alg = def.jet_algorithm();
  def = { alg, std::atof(arg+r) };

  if (prev_alg != undefined_jet_algorithm) {
    std::cout << "\n\033[31mWarning: "
                 "More than 1 argument interpreted as jet algorithm\n"
              << def.description() << "\nreplaced by\n";
    std::cout << def.description() << "\033[0m\n" << std::endl;
  }

  return true;
}

}

#endif
