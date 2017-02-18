#ifndef IVANP_FJJETDEF_FROM_STRING_HH
#define IVANP_FJJETDEF_FROM_STRING_HH

#include <cstdlib>
#include <cctype>
#include <stdexcept>

auto make_jetdef(const char* arg) {
  char c;
  int i = 0;
  for ( ; c = arg[i]; ++i)
    if (c=='\0' || !std::isalpha(c)) break;

  auto cmp = [arg](const char* str, int len) {
    if (i!=len) return false;
    for (int j=0; j<len; ++j)
      if (std::tolower(arg[j]) != str[j]) return false;
    return true;
  };

  using namespace
#ifdef __FASTJET_JETDEFINITION_HH__
  fastjet;
#elif defined(__FJCORE_JETDEFINITION_HH__)
  fjcore;
#else
#error "JetDefinition requires fastjet or fjcore headers included"
#endif

  JetAlgorithm alg;

  if (!std::isdigit(c)) goto x;
  if (cmp("kt",2)) alg = kt_algorithm;
  else if (cmp("antikt",6) || cmp("akt",3)) alg = antikt_algorithm;
  else if (cmp("cambridge",9) || cmp("ca",2)) alg = cambridge_algorithm;
  else goto x;

  const int r = i;
  for ( ; c = arg[i]; ++i) {
    if (c=='\0') break;
    if (!std::isdigit(c)) goto x;
  }

  return JetDefinition(alg,std::atof(arg+r));

x: throw std::invalid_argument(arg);
}

