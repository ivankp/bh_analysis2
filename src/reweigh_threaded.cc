// Written by Ivan Pogrebnyak

#include <iostream>
#include <algorithm>
#include <memory>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include <functional>
#include <future>

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

  static entry current;
};
entry entry::current;

/*
class guarded_pdf {
  LHAPDF::PDF *pdf;
  class wrap {
    static std::mutex guard;
    LHAPDF::PDF *pdf;
  public:
    wrap(LHAPDF::PDF* pdf): pdf(pdf) { guard.lock(); }
    ~wrap() { guard.unlock(); }
    auto operator->() const noexcept { return pdf; }
  };
public:
  ~guarded_pdf() { delete pdf; }
  void operator=(LHAPDF::PDF* pdf) noexcept { this->pdf = pdf; }
  wrap operator->() const noexcept { return { pdf }; }
} pdf;
std::mutex guarded_pdf::wrap::guard;
*/

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
  reweighter(LHAPDF::PDF* pdf)
  : pdf(pdf),
    scale_values(scale_fcns.size()), new_weights(scales.size()),
    _fac_vars(scales_fac.size()), _ren_vars(scales_ren.size())
  { }
  reweighter(const reweighter& o) = delete;
  reweighter(reweighter&& o)
  : pdf(o.pdf),
    scale_values(std::move(o.scale_values)), new_weights(std::move(o.new_weights)),
    _fac_vars(std::move(o._fac_vars)), _ren_vars(std::move(o._ren_vars))
  { o.pdf = nullptr; }

  ~reweighter() { delete pdf; }
  inline void refresh() noexcept {
    memcpy(static_cast<entry*>(this),&entry::current,sizeof(entry));
  }

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

class thread_loop {
  TTree *tree;
  ivanp::timed_counter<Long64_t> &cnt;
  reweighter rew;

  static std::mutex read_mx;

public:
  thread_loop(TTree* tree, decltype(cnt) cnt)
  : tree(tree), cnt(cnt), rew(LHAPDF::mkPDF("CT10nlo",0)) { }

  void operator()() {
    for ( ; ; ) {
      { std::lock_guard<std::mutex> lock(read_mx);
        if (!cnt) break;
        tree->GetEntry(cnt++);
        rew.refresh();
      }
      rew();
    }
  }
  thread_loop(thread_loop&& o)
  : tree(o.tree), cnt(o.cnt), rew(std::move(o.rew)) { };
};
std::mutex thread_loop::read_mx;

double Ht_Higgs(const entry& e) noexcept {
  double _Ht = 0.;
  for (Int_t i=0; i<e.nparticle; ++i) {
    double pt2 = sq(e.px[i],e.py[i]);
    if (e.kf[i]==25) pt2 += sq(125.); // mH^2
    _Ht += std::sqrt(pt2);
  }
  return _Ht;
}


template <typename T>
void branch(TChain& chain, const char* name, T* addr) {
  chain.SetBranchStatus (name, true);
  chain.SetBranchAddress(name, addr);
  chain.AddBranchToCache(name, true);
}

int main(int argc, char* argv[]) {
  if (argc==1) {
    cout << "usage: " << argv[0] << " input.root ..." << endl;
    return 0;
  }

  // pdf = LHAPDF::mkPDF("CT10nlo",0);

  // Open input ntuples root file ===================================
  TChain chain("t3");
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (int i=1; i<argc; ++i) {
    cout << "  " << argv[i] << endl;
    if (!chain.Add(argv[i],0)) return 1;
  }
  chain.SetCacheSize(sq(1024)/4);
  cout << endl;

  chain.SetBranchStatus("*",0);
  branch(chain, "nparticle",    &entry::current.nparticle);
  branch(chain, "px",            entry::current.px);
  branch(chain, "py",            entry::current.py);
  branch(chain, "kf",            entry::current.kf);
  branch(chain, "alphas",       &entry::current.alphas);
  branch(chain, "weight2",      &entry::current.weight2);
  branch(chain, "me_wgt",       &entry::current.me_wgt);
  branch(chain, "me_wgt2",      &entry::current.me_wgt2);
  branch(chain, "x1",           &entry::current.x[0]);
  branch(chain, "x2",           &entry::current.x[1]);
  branch(chain, "x1p",          &entry::current.xp[0]);
  branch(chain, "x2p",          &entry::current.xp[1]);
  branch(chain, "id1",          &entry::current.id[0]);
  branch(chain, "id2",          &entry::current.id[1]);
  branch(chain, "fac_scale",    &entry::current.fac_scale);
  branch(chain, "ren_scale",    &entry::current.ren_scale);
  branch(chain, "usr_wgts",      entry::current.usr_wgts);
  branch(chain, "alphasPower",  &entry::current.alphas_power);
  branch(chain, "part",          entry::current.part);
  chain.StopCacheLearningPhase();

  for (double x : {0.05,0.07,0.10,0.12,0.15,0.20, 0.25,0.3,0.4,0.5})
    scale_fcns.emplace_back([x](const entry& e){ return Ht_Higgs(e)*x; });

  scales_fac = {0,1,2,3,4,5};
  scales_ren = {5,6,7,8,9};
  for (unsigned f=0; f<scales_fac.size(); ++f)
    for (unsigned r=0; r<scales_ren.size(); ++r)
      scales.emplace_back(f,r);

/*
  reweighter rew;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(chain.GetEntries()); ent < 6; ++ent) {
    chain.GetEntry(ent);
    cout << endl;

    rew.refresh();
    rew();

    // test( rew[0] )
  }
*/

  // ivanp::timed_counter<Long64_t> cnt(chain.GetEntries());
  std::vector<std::thread> threads;
  threads.reserve(1);
  ivanp::timed_counter<Long64_t> cnt(10'000);
  for (unsigned i=1; i; --i)
    threads.emplace_back(thread_loop(&chain,cnt));

  // join threads so that program doesn't quit before they are done
  for (auto& thread : threads) thread.join();

  // delete pdf;
  return 0;
}

