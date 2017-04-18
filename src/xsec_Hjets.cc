// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstdio>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "timed_counter.hh"
#include "float_or_double_reader.hh"
#include "parse_args.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

namespace fj = fastjet;

struct ww2 {
  double wtmp = 0., w = 0., w2 = 0.;
  inline void operator+=(double x) noexcept { wtmp += x; }
  void join() noexcept {
    w += wtmp;
    w2 += wtmp*wtmp;
    wtmp = 0.;
  }
};

int main(int argc, char* argv[]) {
  // parse program arguments ========================================
  std::vector<const char*> ntuples, weights;
  const char* tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned need_njets = 0;
  bool w_arg = false;

  for (int i=1; i<argc; ++i) {
    using namespace parse;
    if (!std::strcmp(argv[i],"w")) { w_arg = true; continue; }
    if (!std::strcmp(argv[i],"bh")) { w_arg = false; continue; }
    if (root_file_name(argv[i], w_arg ? weights : ntuples )) continue;
    if (num_jets(argv[i],need_njets)) continue;
    if (jetdef(argv[i],jet_def)) continue;
    if (parse::tree_name(argv[i],tree_name)) continue;

    cerr << "\033[31mUnexpected argument:\033[0m " << argv[i] << endl;
    return 1;
  }
  if (!ntuples.size()) {
    cerr << "\033[31mNo input root files provided\033[0m" << endl;
    return 1;
  }
  if (!need_njets) {
    cerr << "\033[31mNumber of jets is not specified\033[0m" << endl;
    return 1;
  }
  if (jet_def.jet_algorithm() == fj::undefined_jet_algorithm)
    jet_def = fj::JetDefinition(fj::antikt_algorithm,0.4);
  // ================================================================

  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (const char* name : ntuples) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  boost::optional<TChain> weights_chain;
  if (weights.size()) {
    cout << "\n\033[36mInput weights\033[0m:" << endl;
    weights_chain.emplace("weights");
    for (const char* name : weights) {
      if (!weights_chain->Add(name,0)) return 1;
      cout << "  " << name << endl;
    }
    cout << endl;
  }

  // Set up branches for reading
  if (weights_chain) chain.AddFriend(&*weights_chain);
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
  }

  // handle multiple weights
  std::vector<float_or_double_value_reader> _weights;
  if (weights_chain) {
    const TObjArray *bb = weights_chain->GetListOfBranches();
    _weights.reserve(bb->GetEntriesFast()+1);
    _weights.emplace_back(reader,"weight");
    for (const auto* b : *bb)
      _weights.emplace_back(reader,static_cast<const TBranch*>(b)->GetName());
  } else {
    _weights.emplace_back(reader,"weight");
  }
  const auto nw = _weights.size();
  cout << "\n\033[36mWeights\033[0m:\n";
  for (const auto& w : _weights)
    cout << "  " << w.GetBranchName() << '\n';
  cout << endl;

  std::vector<fj::PseudoJet> particles;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n" << endl;

  size_t
    ncount_total = 0,      ncount_selected = 0,
    num_events_total = 0,  num_events_selected = 0,
    num_entries_total = 0, num_entries_selected = 0;

  std::vector<double> this_weight(nw);
  std::vector<ww2> total_weight(nw), selected_weight(nw);

  int current_id = -1, prev_id;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    for (unsigned i=0; i<nw; ++i) this_weight[i] = *_weights[i]; // get weight

    // Keep track of multi-entry events -----------------------------
    prev_id = current_id;
    current_id = *_id;
    ++num_entries_total;
    if (prev_id != current_id) {
      ncount_total += (_ncount ? **_ncount : 1);
      ++num_events_total;
      for (auto& w : total_weight) w.join();
      for (auto& w : selected_weight) w.join();
    }
    // --------------------------------------------------------------

    for (unsigned i=0; i<nw; ++i) total_weight[i] += this_weight[i];

    const size_t np = *_nparticle;
    particles.clear();

    // Read particles -----------------------------------------------
    for (size_t i=0; i<np; ++i) {
      if (_kf[i] != 25) // skip Higgs
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
    }
    // --------------------------------------------------------------

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(jet_pt_cut); // apply pT cut
    // resulting number of jets
    unsigned njets = fj_jets.size();
    // apply eta cut
    for (const auto& jet : fj_jets)
      if (std::abs(jet.eta()) > jet_eta_cut) --njets;
    // --------------------------------------------------------------

    // Apply event cuts ---------------------------------------------
    if (njets < need_njets) continue; // at least needed number of jets
    // --------------------------------------------------------------

    // Keep track of multi-entry events -----------------------------
    ++num_entries_selected;
    if (prev_id != current_id) {
      ncount_selected += (_ncount ? **_ncount : 1);
      ++num_events_selected;
    }
    // --------------------------------------------------------------

    for (auto i=nw; i;) { --i; selected_weight[i] += this_weight[i]; }

  } // END EVENT LOOP
  for (auto& w : total_weight) w.join();
  for (auto& w : selected_weight) w.join();
  // ================================================================

  cout << "\nNumbers" << endl;
  printf("%-11s %-11s  %-11s\n","","total","selected");
  printf("%-11s %11lu  %11lu\n","ncount",ncount_total,ncount_selected);
  printf("%-11s %11lu  %11lu\n","num_events",num_events_total,num_events_selected);
  printf("%-11s %11lu  %11lu\n","num_entries",num_entries_total,num_entries_selected);

  cout << "\nWeights" << endl;
  printf("%-30s %-14s %-14s %-14s %-14s\n",
         "","total","sumw2","selected","sumw2");
  for (unsigned i=0; i<nw; ++i) printf(
    "%-30s %14.8e %14.8e %14.8e %14.8e\n",
    _weights[i].GetBranchName(),
    total_weight[i].w, total_weight[i].w2,
    selected_weight[i].w, selected_weight[i].w2
  );

  cout << "\nCross sections" << endl;
  printf("%-30s %-14s %-14s %-14s %-14s\n",
         "","total","unc","selected","unc");
  for (unsigned i=0; i<nw; ++i) {
    const double xsec = selected_weight[i].w/ncount_total;
    const double unc  = std::sqrt(selected_weight[i].w2)/ncount_total;
    printf(
      "%-30s %14.8e %14.8e %14.8e %14.8e %8.5f%%\n",
      _weights[i].GetBranchName(),
      total_weight[i].w/ncount_total,
      std::sqrt(total_weight[i].w2)/ncount_total,
      xsec, unc, (100.*unc/xsec)
    );
  }
  cout << endl;

  return 0;
}

