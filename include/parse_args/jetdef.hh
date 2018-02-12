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
  unsigned i = 0;
  for ( ; (c = arg[i]); ++i)
    if (c=='\0' || !std::isalpha(c)) break;

  auto cmp = [arg,i](const char* str) {
    const unsigned len = strlen(str);
    if (i!=len) return false;
    for (unsigned j=0; j<len; ++j)
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
  if (cmp("kt")) alg = kt_algorithm;
  else if (cmp("antikt") || cmp("akt")) alg = antikt_algorithm;
  else if (cmp("cambridge") || cmp("ca")) alg = cambridge_algorithm;
  else return false;

  const unsigned r = i;
  for ( ; (c = arg[i]); ++i) {
    if (c=='\0') break;
    if (!std::isdigit(c)) return false;
  }

  const JetAlgorithm prev_alg = def.jet_algorithm();
  def = { alg, std::atof(arg+r)*0.1 };

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
