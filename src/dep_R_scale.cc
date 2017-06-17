#include <iostream>
#include <fstream>

#include <TChain.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "reweighter.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;

namespace fj = fastjet;

using namespace ivanp::math;
#include "scales.hh"

struct ww2 {
  double wtmp = 0., w = 0., w2 = 0.;
  inline void operator+=(double x) noexcept { wtmp += x; }
  void join() noexcept {
    w += wtmp;
    w2 += wtmp*wtmp;
    wtmp = 0.;
  }
};
std::ostream& operator<<(std::ostream &os, const ww2& wc) {
  os << setw(16) << wc.w << setw(15) << wc.w2 << '\n';
  return os;
}

// int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " out.root in.root ..." << endl;
    return 1;
  }

  // Open input ntuples root file ===================================
  TChain chain("t3");
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (const char* f : ntuples) {
    cout << "  " << f << endl;
    if (!chain.Add(f,0)) return 1;
  }

  scale_defs sd;

  std::vector<double> Ht_fracs {1,2,3,4,6,8,12,16};
  for (double x : Ht_fracs) Ht_fracs.push_back(1./x);
  std::sort(Ht_fracs.begin(),Ht_fracs.end());

  for (double x : Ht_fracs) {
    static unsigned i = 0;
    sd.scale_fcns.emplace_back([x](const entry& e){ return HT_hat_pp(e)*x; });
    sd.scales_fac.emplace_back(i);
    sd.scales_ren.emplace_back(i);
    sd.scales.emplace_back(i,i);
    ++i;
  }

  reweighter rew(chain,sd,"CT10nlo");

  Int_t event_id, prev_id = -1;
  branch(chain, "id", &event_id);

  std::vector<ww2> scales_weights(sd.scales.size());

  const entry& e = *reinterpret_cast<const entry*>(&rew);
  std::vector<fj::PseudoJet> particles;

  unsigned n_events_total = 0, n_events_selected = 0;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n";
  test(pt_cut)
  cout << endl;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(chain.GetEntries()); !!ent; ++ent) {
    chain.GetEntry(ent);

    const bool new_event = (event_id != prev_id);
    if (new_event) {
      ++n_events_total;
      prev_id = event_id;
      for (auto& w : scales_weights) w.join();
    }

    // Read particles -----------------------------------------------
    particles.clear();
    for (int i=0; i<e.nparticle; ++i) {
      if (e.kf[i] != 25) // skip Higgs
        particles.emplace_back(e.px[i],e.py[i],pz[i],E[i]);
    }
    // --------------------------------------------------------------

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(pt_cut); // apply pT cut
    // resulting number of jets
    unsigned njets = fj_jets.size();
    // apply eta cut
    for (const auto& jet : fj_jets)
      if (std::abs(jet.eta()) > 4.4) --njets;
    // --------------------------------------------------------------

    // Apply event cuts ---------------------------------------------
    if (njets < need_njets) continue; // at least needed number of jets
    if (new_event) ++n_events_selected;
    // --------------------------------------------------------------

    rew();

    for (auto i=sd.scales.size(); i;) { --i; scales_weights[i] += rew[i]; }
  }
  for (auto& w : scales_weights) w.join();

  std::ofstream fout(dat);
  if (!fout) {
    cerr << "cannot open output file " << dat << endl;
    return 1;
  }

  fout << std::scientific << std::setprecision(8);
  fout << n_events_total <<' '<< n_events_selected << '\n';
  for (const auto& wc : scales_weights) {
    static int i = 0;
    fout << setw(14) << Ht_fracs[sd.scales_fac[*sd.scales[i].fac]];
    fout << setw(15) << Ht_fracs[sd.scales_ren[*sd.scales[i].ren]];
    fout << wc;
    ++i;
  }

  return 0;
}

