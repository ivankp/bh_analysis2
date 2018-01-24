#ifndef MULTIBIN_HH
#define MULTIBIN_HH

#include "multibin.hh"

struct multibin {
  static std::vector<double> weights;
  static unsigned wi;
};
std::vector<double> multibin::weights;
unsigned multibin::wi;

#endif
