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
#include <TDirectory.h>
#include <TKey.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "binner.hh"
#include "re_axes.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"
#include "slice.hh"

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
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

template <typename A, typename B>
TH1D* make_TH1D(const char* name, const A& axis, B begin, B end) {
  const auto edges = ivanp::make_vector_of_edges<double>(axis);
  TH1D* h = new TH1D(name,"",edges.size()-1,edges.data());
  h->Sumw2();
  TArrayD& sumw2 = *h->GetSumw2();
  size_t n_total = 0, i = 0;
  for (auto bin=begin; bin!=end; ++bin, ++i) {
    h->SetBinContent(i,bin->w);
    sumw2[i] = bin->w2;
    n_total += bin->n;
  }
  h->SetEntries(n_total);
  return h;
}

template <typename A>
TH1D* root_hist(const ivanp::binner<hist_bin,std::tuple<A>>& h,
  const std::string& name
) {
  return make_TH1D(name.c_str(),h.axis(),h.bins().begin(),h.bins().end());
}

int main(int argc, char* argv[]) {
  const fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
  constexpr unsigned need_njets = 2;

  // Define histograms ==============================================
  size_t ncount = 0, num_events = 0, num_selected = 0;

  hist<int> h_Njets({need_njets+2,0,need_njets+2});

  re_axes ra("example.bins");
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
  TTreeReaderValue<Double_t> _weight(reader,"weight");

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
  TTreeReaderValue<Int_t> _id1(reader,"id1");
  TTreeReaderValue<Int_t> _id2(reader,"id2");

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1, curr_id;

  fastjet::ClusterSequence::print_banner(); // get it out of the way

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next() /* && ent < 1 */; ++ent) {
    hist_bin::weight = *_weight; // Read weight

    // Keep track of multi-entry events -----------------------------
    curr_id = *_id;
    if (prev_id!=curr_id) {
      prev_id = curr_id;
      ncount += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------

    particles.clear();
    particles.reserve(*_nparticle);
    boost::optional<TLorentzVector> higgs;

    // Read particles -----------------------------------------------
    for (size_t i=0, n=_kf.GetSize(); i<n; ++i) {
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
      if (abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
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
  // ****************************************************************

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(argv[1],"recreate");
  if (fout->IsZombie()) return 1;

  fout->mkdir("weight_JetAntiKt4")->cd();

  // write root historgrams
  root_hist(h_Njets,"jets_N_excl");
  h_Njets.integrate_left();
  root_hist(h_Njets,"jets_N_incl");

  for (auto& h : re_hist<1>::all) root_hist(*h,h.name);

  fout->cd();
  (new TH1D("N","N",1,0,1))->SetBinContent(1,ncount);
  fout->Write();
  cout << "\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}

