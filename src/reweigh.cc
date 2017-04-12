// Written by Ivan Pogrebnyak

#include <iostream>
#include <algorithm>
#include <memory>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <functional>

#include <boost/optional.hpp>

#include <TChain.h>

#include "LHAPDF/LHAPDF.h"

#include "timed_counter.hh"
// #include "catstr.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

// using ivanp::cat;

constexpr std::array<int,10> quarks {
   1,-1, 2,-2, 3,-3, 4,-4, 5,-5
};

struct entry {
  Int_t           nparticle;
  Float_t         px[8]; //[nparticle]
  Float_t         py[8]; //[nparticle]
  // Float_t         pz[8]; //[nparticle]
  // Float_t         E [8]; //[nparticle]
  Double_t        alphas;
  Double_t        weight;
  Double_t        weight2;
  Double_t        me_wgt;
  Double_t        me_wgt2;
  Double_t        x[2];
  Double_t        xp[2];
  Int_t           id[2];
  Double_t        fac_scale;
  Double_t        ren_scale;
  Double_t        usr_wgts[18]; //[nuwgt]
  Short_t         alphas_power;
  Char_t          part[2];

  static entry current;
};
entry entry::current;

LHAPDF::PDF *pdf;

std::vector<std::function<double(const entry&)>> scale_fcns;
std::vector<unsigned> scales_fac, scales_ren;
struct scale_def {
  boost::optional<unsigned> fac, ren;
  scale_def(unsigned f, unsigned r): fac(f), ren(r) { }
};
std::vector<scale_def> scales;

class reweigh : entry { // thread instance for an entry
  std::vector<double> scale_values, new_weights;

  struct fac_vars { double m, ff; };
  struct ren_vars { double ar, w0; };
  std::vector<fac_vars> _fac_vars;
  std::vector<ren_vars> _ren_vars;

  double fr1(unsigned r, double muF) const {
    const double x = this->x[r];
    if (id[r] != 21) return pdf->xfxQ(id[r], x, muF)/x;
    else {
      double f = 0.;
      for (int q : quarks) f += pdf->xfxQ(q, x, muF)/x;
      return f;
    }
  };

  double fr2(unsigned r, double muF) const {
    const double x  = this->x[r];
    const double xp = x/this->xp[r];
    if (id[r] != 21) return pdf->xfxQ(id[r], xp, muF)/x;
    else {
      double f = 0.;
      for (int q : quarks) f += pdf->xfxQ(q, xp, muF)/x;
      return f;
    }
  };

  double fr3(unsigned r, double muF) const {
    return pdf->xfxQ(21, x[r], muF)/x[r];
  }
  double fr4(unsigned r, double muF) const {
    return pdf->xfxQ(21, x[r]/xp[r], muF)/x[r];
  }

  void fac_calc(unsigned i) {
    fac_vars& v = _fac_vars[i];
    const double muF = scale_values[scales_fac[i]];

    double f[2];
    for (int i=0; i<2; ++i) f[i] = pdf->xfxQ(id[i], x[i], muF)/x[i];
    v.ff = f[0]*f[1];

    if (part[0]=='I') {
      const double lf = 2.*std::log(muF/fac_scale);
      double w[8];
      for (int i=0; i<8; ++i) w[i] = usr_wgts[i+2] + usr_wgts[i+10]*lf;

      v.m = ( fr1(0,muF)*w[0] + fr2(0,muF)*w[1]
            + fr3(0,muF)*w[2] + fr4(0,muF)*w[3] )*f[1]
          + ( fr1(1,muF)*w[4] + fr2(1,muF)*w[5]
            + fr3(1,muF)*w[6] + fr4(1,muF)*w[7] )*f[0];
    } else v.m = 0.;
  }

  void ren_calc(unsigned i) {
    ren_vars& v = _ren_vars[i];
    const double muR = scale_values[scales_ren[i]];

    v.ar = std::pow(pdf->alphasQ(muR)/alphas, alphas_power);

    if (part[0]=='V' || part[0]=='I') {
      const double lr = 2.*std::log(muR/ren_scale);
      v.w0 = me_wgt + lr*usr_wgts[0] + 0.5*lr*lr*usr_wgts[0];
    } else {
      v.w0 = me_wgt2;
    }
  }

public:
  reweigh()
  : scale_values(scale_fcns.size()), new_weights(scales.size()),
    _fac_vars(scales_fac.size()), _ren_vars(scales_ren.size())
  { }
  inline void refresh() noexcept {
    memcpy(static_cast<entry*>(this),&entry::current,sizeof(entry));
  }

  void operator()() {
    for (unsigned i=0; i<scale_fcns.size(); ++i)
      scale_values[i] = scale_fcns[i](*this);

    for (auto i : scales_fac) fac_calc(i);
    for (auto i : scales_ren) ren_calc(i);

    for (unsigned i=0; i<scales.size(); ++i) { // combine
      const scale_def& scale = scales[i];
      const fac_vars& fac = _fac_vars[i];
      const ren_vars& ren = _ren_vars[i];

      const double m = scale.fac ? (fac.m + ren.w0*fac.ff) : weight2;

      new_weights[i] = scale.ren ? ren.ar * m : m;
    }
  }
  inline double operator[](unsigned i) const { return new_weights[i]; }
};

template <typename T>
void branch(TChain& chain, const char* name, T* addr) {
  chain.SetBranchStatus (name, true);
  chain.SetBranchAddress(name, addr);
}

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " input.root output.root" << endl;
    return 0;
  }

  pdf = LHAPDF::mkPDF("CT10nlo",0);

  // Open input ntuples root file ===================================
  TChain chain("BHSntuples");
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (int i=1; i<argc; ++i) {
    cout << "  " << argv[i] << endl;
    if (!chain.Add(argv[i],0)) return 1;
  }
  cout << endl;

  chain.SetBranchStatus("*",0);
  branch(chain, "alphas",      &entry::current.alphas);
  branch(chain, "weight",      &entry::current.weight);
  branch(chain, "me_wgt",      &entry::current.me_wgt);
  branch(chain, "x1",          &entry::current.x[0]);
  branch(chain, "x2",          &entry::current.x[1]);
  branch(chain, "x1p",         &entry::current.xp[0]);
  branch(chain, "x2p",         &entry::current.xp[1]);
  branch(chain, "id1",         &entry::current.id[0]);
  branch(chain, "id2",         &entry::current.id[1]);
  branch(chain, "fac_scale",   &entry::current.fac_scale);
  branch(chain, "ren_scale",   &entry::current.ren_scale);
  branch(chain, "usr_wgts",     entry::current.usr_wgts);
  branch(chain, "alphasPower", &entry::current.alphas_power);
  branch(chain, "part",         entry::current.part);

  scale_fcns.emplace_back([](const entry&){ return 66.3463; });

  scales_fac.emplace_back(0);
  scales_ren.emplace_back(0);
  scales.emplace_back(0,0);

  reweigh rew;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(chain.GetEntries()); ent < 6; ++ent) {
    chain.GetEntry(ent);
    cout << endl;

    rew.refresh();
    rew();

    test( rew[0] )
  }

  delete pdf;
  return 0;
}

