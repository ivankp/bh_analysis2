// Written by Ivan Pogrebnyak

#include <iostream>
// #include <stdexcept>

#include <TChain.h>

#include "math.hh"
#include "timed_counter.hh"
#include "reweighter.hh"
// #include "catstr.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

// using ivanp::cat;
using namespace ivanp::math;

double Ht_Higgs(const entry& e) noexcept {
  double _Ht = 0.;
  for (Int_t i=0; i<e.nparticle; ++i) {
    // double pt2 = sq(e.px[i],e.py[i]);
    // if (e.kf[i]==25) pt2 += sq(125.); // mH^2
    // _Ht += std::sqrt(pt2);
    _Ht += std::sqrt( e.kf[i]==25
        ? sq(e.E[i])-sq(e.pz[i]) // mT for Higgs
        : sq(e.px[i],e.py[i]) ); // pT
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

  scale_defs sd;

  // for (double x : {0.05,0.07,0.10,0.12,0.15,0.20, 0.25,0.3,0.4,0.5})
  for (double x : {0.5,1.,0.25})
    sd.scale_fcns.emplace_back([x](const entry& e){ return Ht_Higgs(e)*x; });

  // sd.scales_fac = {0,1,2,3,4,5};
  // sd.scales_ren = {5,6,7,8,9};
  sd.scales_fac = {0,1,2};
  sd.scales_ren = {0,1,2};
  // for (unsigned f=0; f<sd.scales_fac.size(); ++f)
  //   for (unsigned r=0; r<sd.scales_ren.size(); ++r)
  //     sd.scales.emplace_back(f,r);
  sd.scales.emplace_back(0,0);
  sd.scales.emplace_back(0,1);
  sd.scales.emplace_back(1,0);
  sd.scales.emplace_back(1,1);
  sd.scales.emplace_back(0,2);
  sd.scales.emplace_back(2,0);
  sd.scales.emplace_back(2,2);

  reweighter rew(chain,(1<<19),"CT10nlo",sd);

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  // for (tc ent(chain.GetEntries()); !!ent; ++ent) {
  for (tc ent(1); !!ent; ++ent) {
    chain.GetEntry(ent);
    // cout << endl;

    rew();

    // test( rew[0] )
    for (unsigned i=0, n=sd.scales.size(); i<n; ++i)
      cout << rew[i] << endl;
  }

  return 0;
}

