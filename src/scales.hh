#ifndef SCALES_DEFINITIONS_HH
#define SCALES_DEFINITIONS_HH

#include <map>

double HT(const entry& e) noexcept {
  double _HT = 0.;
  for (Int_t i=0; i<e.nparticle; ++i)
    _HT += std::sqrt( sq(e.px[i],e.py[i]) ); // pT
  return _HT;
}

double HT_hat_p(const entry& e) noexcept {
  double _HT = 0.;
  for (Int_t i=0; i<e.nparticle; ++i) {
    _HT += std::sqrt( e.kf[i]==25
        ? sq(e.E[i])-sq(e.pz[i]) // mT for Higgs
        : sq(e.px[i],e.py[i]) ); // pT
  }
  return _HT;
}
double HT_hat_p2(const entry& e) noexcept {
  return HT_hat_p(e)*0.5;
}

double HT_hat_pp(const entry& e) noexcept {
  return 125. + 0.5*HT(e);
}

const std::map<std::string,double(*)(const entry& e)> scales {
  {"HT", &HT},
  {"HT1", &HT_hat_p2},
  {"HT2", &HT_hat_pp}
};
const std::map<std::string,std::string> scales_pretty {
  {"HT", "HT"},
  {"HT1", "ĤT\'/2"},
  {"HT2", "ĤT\'\'"}
};

#endif
