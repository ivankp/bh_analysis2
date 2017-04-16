#include <vector>
#include <functional>

#include <boost/optional.hpp>

#include <TTree.h>

template <typename T>
void branch(TTree& tree, const char* name, T* addr) {
  tree.SetBranchStatus (name, true);
  tree.SetBranchAddress(name, addr);
  tree.AddBranchToCache(name, true);
}

struct entry {
  Int_t           nparticle;
  Float_t         px[8]; //[nparticle]
  Float_t         py[8]; //[nparticle]
  Int_t           kf[8];
  Double_t        alphas;
  Double_t        weight2;
  Double_t        me_wgt;
  Double_t        me_wgt2;
  Double_t        x[2];
  Double_t        xp[2];
  Int_t           id[2];
  Double_t        fac_scale;
  Double_t        ren_scale;
  Double_t        usr_wgts[18]; //[nuwgt]
  Char_t          alphas_power;
  Char_t          part[2];

  entry(TTree& tree, Long64_t cacheSize);
};

struct scale_defs {
  std::vector<std::function<double(const entry&)>> scale_fcns;
  std::vector<unsigned> scales_fac, scales_ren;
  struct scale {
    boost::optional<unsigned> fac, ren;
    scale(unsigned f, unsigned r): fac(f), ren(r) { }
    scale(boost::none_t f, unsigned r): fac(f), ren(r) { }
    scale(unsigned f, boost::none_t r): fac(f), ren(r) { }
  };
  std::vector<scale> scales;
};

namespace LHAPDF {
  class PDF;
}

class reweighter : entry {
  LHAPDF::PDF *pdf; // owned here
  const scale_defs& sd;

  std::vector<double> scale_values, new_weights;

  struct fac_vars { double m, ff; };
  struct ren_vars { double ar, w0; };
  std::vector<fac_vars> _fac_vars;
  std::vector<ren_vars> _ren_vars;

public:
  reweighter(TTree& tree, Long64_t cacheSize,
             const char* pdf_name, const scale_defs& sd);
  ~reweighter();

  double fr1(unsigned r, double muF) const;
  double fr2(unsigned r, double muF) const;
  double fr3(unsigned r, double muF) const;
  double fr4(unsigned r, double muF) const;

  void fac_calc(unsigned i);
  void ren_calc(unsigned i);

public:
  void operator()();
  inline double operator[](unsigned i) const { return new_weights[i]; }
};

