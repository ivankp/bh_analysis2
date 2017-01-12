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
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "binner.hh"
#include "timed_counter.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using std::cout;
using std::cerr;
using std::endl;

using p4_t = Float_t;

namespace fj = fastjet;

struct hist_bin {
  double w, w2;
  size_t n;
  hist_bin(): w(0.), w2(0.), n(0) { }
  inline void operator+=(double weight) {
    w += weight;
    w2 += weight*weight;
    ++n;
  }
};

struct ntuple_filler {
  static double weight;
  template <typename T>
  void operator()(T& bin) noexcept { bin += weight; }
};
double ntuple_filler::weight;

template <typename T>
using hist = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<ivanp::uniform_axis<T>>>,
  std::vector<hist_bin>,
  ntuple_filler>;

template <typename T>
TH1D* root_hist(const hist<T>& h, const char* name) {
  TH1D* hr = new TH1D(name,"",h.axis().nbins(),h.axis().min(),h.axis().max());
  hr->Sumw2();
  TArrayD& sumw2 = *hr->GetSumw2();
  const auto& bins = h.bins();
  size_t n_total = 0;
  for (int i=0, n=bins.size(); i<n; ++i) {
    const auto& bin = bins[i];
    hr->SetBinContent(i,bin.w);
    sumw2[i] = bin.w2;
    n_total += bin.n;
  }
  hr->SetEntries(n_total);
  return hr;
}

int main(int argc, char* argv[])
{
  fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
  const unsigned need_njets = 2;

  // Define histograms ==============================================
  size_t N = 0, num_events = 0, num_selected = 0;

  hist<int> h_Njets_incl({4,0,4}), h_Njets_excl({4,0,4});
  hist<double> h_H_pT({70,0,700}), h_H_y({36,-4.5,4.5}), h_H_eta({36,-4.5,4.5}),
  h_H_phi({36,-M_PI,M_PI}), h_H_mass({50,105,155});
  // ================================================================

  // Open input ntuple root file
  auto fin = std::make_unique<TFile>(argv[1],"read");
  if (fin->IsZombie()) return 1;

  // Set up branches for reading
  TTreeReader reader("t3",fin.get());
  if (!reader.GetTree()) {
    cerr << "\033[31mNo tree \"t3\" in file"
         << argv[1] <<"\033[0m"<< endl;
    return 1;
  }

  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");
  TTreeReaderArray<p4_t>  _px(reader,"px");
  TTreeReaderArray<p4_t>  _py(reader,"py");
  TTreeReaderArray<p4_t>  _pz(reader,"pz");
  TTreeReaderArray<p4_t>  _E (reader,"E");
  TTreeReaderValue<Double_t> _weight(reader,"weight");
  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  // boost::optional<TTreeReaderValue<Char_t>> _part;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
    // else if (!strcmp(bo->GetName(),"part")) _part.emplace(reader,"part");
  }

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1, curr_id;

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {

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

    curr_id = *_id;
    if (prev_id!=curr_id) {
      prev_id = curr_id;
      N += ( _ncount ? **_ncount : 1);
      ++num_events;
    }

    // Cluster jets -------------------------------------------------
    auto jets = fj::ClusterSequence(particles,jet_def)
      .inclusive_jets(jet_pt_cut);
    // if (jets.size() < need_njets) continue;
    std::remove_if( jets.begin(), jets.end(),
      [=](const fj::PseudoJet& j){ return (j.eta() > jet_eta_cut); });
    // if (jets.size() < need_njets) continue;
    std::sort( jets.begin(), jets.end(),
      [](const fj::PseudoJet& a, const fj::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    const unsigned njets = jets.size();
    // --------------------------------------------------------------

    // Define variables ---------------------------------------------
    ntuple_filler::weight = *_weight;

    const double H_pT = higgs->Pt();
    // --------------------------------------------------------------

    ++num_selected;

    // --------------------------------------------------------------
    for (unsigned j=1; j<=njets+1; ++j) h_Njets_incl.fill_bin(j);
    h_Njets_excl.fill_bin(njets+1); // nj+1 because nj==0 is bin 1

    h_H_pT(H_pT);
    h_H_y(higgs->Rapidity());
    h_H_eta(higgs->Eta());
    h_H_phi(higgs->Phi());
    h_H_mass(higgs->M());

    // --------------------------------------------------------------
  } // END EVENT LOOP
  // ****************************************************************

  test(N)

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;

  // Open input ntuple root file
  auto fout = std::make_unique<TFile>(argv[2],"recreate");
  if (fout->IsZombie()) return 1;

  // write root historgrams
#define rh(name) root_hist(h_##name,#name);
  rh(Njets_incl)
  rh(Njets_excl)

  rh(H_pT)
  rh(H_y)
  rh(H_eta)
  rh(H_phi)
  rh(H_mass)

  fout->Write();

  return 0;
}

