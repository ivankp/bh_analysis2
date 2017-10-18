#ifndef MULTIBIN_HH
#define MULTIBIN_HH

#include "multibin.hh"

struct multibin {
  static std::vector<double> weights;
  static unsigned wi;
  static int current_id;
};
std::vector<double> multibin::weights;
unsigned multibin::wi;
int multibin::current_id;

#endif
