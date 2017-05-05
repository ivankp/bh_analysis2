#ifndef IVANP_HIGGS2DIPHOTON_HH
#define IVANP_HIGGS2DIPHOTON_HH

#include <utility>
#include <random>

#include <TLorentzVector.h>

class Higgs2diphoton {
  std::mt19937 gen; // mersenne twister random number generator
  std::uniform_real_distribution<double> phi_dist; // φ
  std::uniform_real_distribution<double> cts_dist; // cos(θ*)

public:
  Higgs2diphoton();

  std::pair<TLorentzVector,TLorentzVector>
  operator()(const TLorentzVector& Higgs);
};


#endif
