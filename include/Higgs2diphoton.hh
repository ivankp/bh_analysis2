#ifndef IVANP_HIGGS2DIPHOTON_HH
#define IVANP_HIGGS2DIPHOTON_HH

#include <utility>
#include <random>

#include <TLorentzVector.h>

class Higgs2diphoton {
  std::mt19937 gen; // mersenne twister random number generator
  std::uniform_real_distribution<double> phi_dist; // φ
  std::uniform_real_distribution<double> cts_dist; // cos(θ*)

  double phi, cts, sts, cos_phi, sin_phi; // remember last values

public:
  Higgs2diphoton();

  using seed_type = typename decltype(gen)::result_type;
  Higgs2diphoton(seed_type seed);

  std::pair<TLorentzVector,TLorentzVector>
  operator()(const TLorentzVector& Higgs, bool new_kin=true);
};


#endif
