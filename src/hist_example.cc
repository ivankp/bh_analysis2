// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstring>
#include <algorithm>
#include <memory>

#include <boost/optional.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "re_axes.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;
using ivanp::reserve;

namespace fj = fastjet;
using namespace ivanp::math;

struct Jet {
  TLorentzVector p;
  double pT, y, eta, phi, mass;
  template <typename P>
  Jet(const P& _p): p(_p[0],_p[1],_p[2],_p[3]),
    pT(p.Pt()), y(p.Rapidity()), eta(p.Eta()), phi(p.Phi()), mass(p.M())
  { }
};

struct hist_bin {
  static double weight;
  double w, w2;
  size_t n;
  hist_bin(): w(0.), w2(0.), n(0) { }
  inline void operator++() noexcept {
    w += weight;
    w2 += weight*weight;
    ++n;
  }
  inline hist_bin& operator+=(const hist_bin& b) noexcept {
    w += b.w;
    w2 += b.w2;
    n += b.n;
    return *this;
  }
};
double hist_bin::weight;

template <typename... Axes>
using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <bool... OF>
using re_hist = ivanp::binner<hist_bin, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;

int main(int argc, char* argv[]) {
  const fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
  const unsigned need_njets = 2;

  // Define histograms ==============================================
  hist<int> h_Njets({need_njets+2,0,need_njets+2});

  re_axes ra("config/example.bins");
#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  a_(y) a_(phi)

  h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

  auto h_jet_pT   = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_y    = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_eta  = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_phi  = reserve<re_hist<1>>(need_njets+1);
  auto h_jet_mass = reserve<re_hist<1>>(need_njets+1);

  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_pT");
    h_jet_pT.emplace_back(name,ra[name]);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_y");
    h_jet_y.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_eta");
    h_jet_eta.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_phi");
    h_jet_phi.emplace_back(name,a_phi);
  }
  for (unsigned i=0; i<need_njets+1; ++i) {
    auto name = cat("jet",i+1,"_mass");
    h_jet_mass.emplace_back(name,ra[name]);
  }

  // ================================================================

  // Open input ntuple root file
  TChain chain("t3");
  for (int i=2; i<argc; ++i) {
    if (!chain.Add(argv[i],0)) return 1;
    cout << "\033[36mInput\033[0m: " << argv[i] << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");
  TTreeReaderValue<Double_t> _weight(reader,"weight2");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  // boost::optional<TTreeReaderValue<Char_t>> _part;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
    // else if (!strcmp(bo->GetName(),"part")) _part.emplace(reader,"part");
  }

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1, curr_id;

  fastjet::ClusterSequence::print_banner(); // get it out of the way

  size_t ncount = 0, num_events = 0, num_selected = 0;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    hist_bin::weight = *_weight; // Read weight

    // Keep track of multi-entry events -----------------------------
    curr_id = *_id;
    if (prev_id!=curr_id) {
      prev_id = curr_id;
      ncount += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------

    const size_t np = *_nparticle;
    particles.clear();
    particles.reserve(np);
    boost::optional<TLorentzVector> higgs;

    // Read particles -----------------------------------------------
    for (size_t i=0; i<np; ++i) {
      if (_kf[i] == 25) {
        higgs.emplace(_px[i],_py[i],_pz[i],_E[i]);
      } else {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }
    if (!higgs) cerr << "\033[31mNo Higgs in entry " << ent <<"\033[0m"<< endl;
    // --------------------------------------------------------------

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(jet_pt_cut); // apply pT cut
    // apply eta cut
    for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
      if (std::abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
      else ++it;
    }
    // sort by pT
    std::sort( fj_jets.begin(), fj_jets.end(),
      [](const fj::PseudoJet& a, const fj::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    // resulting number of jets
    const unsigned njets = fj_jets.size();
    // --------------------------------------------------------------

    // Fill Njets histograms before cuts
    h_Njets.fill_bin(njets+1); // njets+1 because njets==0 is bin 1
    // --------------------------------------------------------------

    // Apply event cuts ---------------------------------------------
    if (njets < need_njets) continue; // at least needed number of jets
    // --------------------------------------------------------------

    // Define variables ---------------------------------------------
    const double H_pT = higgs->Pt();

    std::vector<Jet> jets; // compute jet variables
    jets.reserve(njets);
    for (const auto& jet : fj_jets) jets.emplace_back(jet);

    double HT = H_pT;
    for (const auto& j : jets) HT += j.pT;

    double H_y = higgs->Rapidity();
    double H_phi = higgs->Phi();
    // --------------------------------------------------------------

    ++num_selected;

    // Fill histograms ----------------------------------------------
    h_HT(HT);

    h_H_pT(H_pT);
    h_H_y(H_y);
    h_H_eta(higgs->Eta());
    h_H_phi(H_phi);
    h_H_mass(higgs->M());

    // jet histograms
    for (unsigned i=0, n=std::min(njets,need_njets+1); i<n; ++i) {
      h_jet_pT  [i](jets[i].pT  );
      h_jet_y   [i](jets[i].y   );
      h_jet_eta [i](jets[i].eta );
      h_jet_phi [i](jets[i].phi );
      h_jet_mass[i](jets[i].mass);
    }
    // --------------------------------------------------------------

  } // END EVENT LOOP
  // ================================================================

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << endl;

  // Open output root file for histograms
  TFile fout(argv[1],"recreate");
  if (fout.IsZombie()) return 1;

  fout.mkdir("weight_JetAntiKt4")->cd();

  // write root historgrams
  using ivanp::root::to_root;
  using ivanp::root::slice_to_root;

  auto* h_Njets_excl = to_root(h_Njets,"jets_N_excl");
  h_Njets.integrate_left();
  auto* h_Njets_incl = to_root(h_Njets,"jets_N_incl");
  h_Njets_incl->SetEntries( h_Njets_excl->GetEntries() );

  for (auto& h : re_hist<1>::all) to_root(*h,h.name);

  fout.cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_events);
  fout.Write();
  cout << "\033[32mOutput\033[0m: " << fout.GetName() << endl;

  return 0;
}

