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

#include "math.hh"
#include "timed_counter.hh"
// #include "catstr.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

// using ivanp::cat;
using namespace ivanp::math;

constexpr std::array<int,10> quarks {
  1,-1, 2,-2, 3,-3, 4,-4, 5,-5
};

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

  entry(TTree& tree, Long64_t cacheSize) {
    tree.SetCacheSize(cacheSize);
    tree.SetBranchStatus("*",0);
    branch(tree, "nparticle",   &nparticle);
    branch(tree, "px",           px);
    branch(tree, "py",           py);
    branch(tree, "kf",           kf);
    branch(tree, "alphas",      &alphas);
    branch(tree, "weight2",     &weight2);
    branch(tree, "me_wgt",      &me_wgt);
    branch(tree, "me_wgt2",     &me_wgt2);
    branch(tree, "x1",          &x[0]);
    branch(tree, "x2",          &x[1]);
    branch(tree, "x1p",         &xp[0]);
    branch(tree, "x2p",         &xp[1]);
    branch(tree, "id1",         &id[0]);
    branch(tree, "id2",         &id[1]);
    branch(tree, "fac_scale",   &fac_scale);
    branch(tree, "ren_scale",   &ren_scale);
    branch(tree, "usr_wgts",     usr_wgts);
    branch(tree, "alphasPower", &alphas_power);
    branch(tree, "part",         part);
    tree.StopCacheLearningPhase();
  }
};

std::vector<std::function<double(const entry&)>> scale_fcns;
std::vector<unsigned> scales_fac, scales_ren;
struct scale_def {
  boost::optional<unsigned> fac, ren;
  scale_def(unsigned f, unsigned r): fac(f), ren(r) { }
  scale_def(boost::none_t f, unsigned r): fac(f), ren(r) { }
  scale_def(unsigned f, boost::none_t r): fac(f), ren(r) { }
};
std::vector<scale_def> scales;

class reweighter : entry { // thread instance for an entry
  LHAPDF::PDF *pdf; // owned here
  std::vector<double> scale_values, new_weights;

  struct fac_vars { double m, ff; };
  struct ren_vars { double ar, w0; };
  std::vector<fac_vars> _fac_vars;
  std::vector<ren_vars> _ren_vars;

public:
  reweighter(TTree& tree, Long64_t cacheSize, LHAPDF::PDF* pdf)
  : entry(tree,cacheSize), pdf(pdf),
    scale_values(scale_fcns.size()), new_weights(scales.size()),
    _fac_vars(scales_fac.size()), _ren_vars(scales_ren.size())
  { }
  ~reweighter() { delete pdf; }

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
      v.w0 = me_wgt + lr*usr_wgts[0] + 0.5*lr*lr*usr_wgts[1];
    } else {
      v.w0 = me_wgt2;
    }
  }

public:
  void operator()() {
    for (unsigned i=0; i<scale_fcns.size(); ++i)
      scale_values[i] = scale_fcns[i](*this);

    for (unsigned i=0; i<scales_fac.size(); ++i) fac_calc(i);
    for (unsigned i=0; i<scales_ren.size(); ++i) ren_calc(i);

    for (unsigned i=0; i<scales.size(); ++i) { // combine
      const scale_def& scale = scales[i];

      double& w = new_weights[i];
      if (scale.fac) {
        const fac_vars& fac = _fac_vars[*scale.fac];
        w = fac.m;
        if (scale.ren) {
          const ren_vars& ren = _ren_vars[*scale.ren];
          w += ren.w0 * fac.ff;
        } else {
          w += me_wgt2 * fac.ff;
        }
      } else {
        w = weight2;
      }
      if (scale.ren) {
        const ren_vars& ren = _ren_vars[*scale.ren];
        w *= ren.ar;
      }

      // test( w )
    }
    // cout << endl;
  }
  inline double operator[](unsigned i) const { return new_weights[i]; }
};

double Ht_Higgs(const entry& e) noexcept {
  double _Ht = 0.;
  for (Int_t i=0; i<e.nparticle; ++i) {
    double pt2 = sq(e.px[i],e.py[i]);
    if (e.kf[i]==25) pt2 += sq(125.); // mH^2
    _Ht += std::sqrt(pt2);
  }
  return _Ht;
}

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " input.root ..." << endl;
    return 0;
  }

  // Open input ntuples root file ===================================
  TChain chain("t3");
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (int i=1; i<argc; ++i) {
    cout << "  " << argv[i] << endl;
    if (!chain.Add(argv[i],0)) return 1;
  }

  for (double x : {0.05,0.07,0.10,0.12,0.15,0.20, 0.25,0.3,0.4,0.5})
    scale_fcns.emplace_back([x](const entry& e){ return Ht_Higgs(e)*x; });

  scales_fac = {0,1,2,3,4,5};
  scales_ren = {5,6,7,8,9};
  for (unsigned f=0; f<scales_fac.size(); ++f)
    for (unsigned r=0; r<scales_ren.size(); ++r)
      scales.emplace_back(f,r);

  reweighter rew(chain,(1<<19),LHAPDF::mkPDF("CT10nlo",0));

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  // for (tc ent(chain.GetEntries()); !!ent; ++ent) {
  for (tc ent(100'000); !!ent; ++ent) {
    chain.GetEntry(ent);
    // cout << endl;

    rew();

    // test( rew[0] )
  }

  return 0;
}

