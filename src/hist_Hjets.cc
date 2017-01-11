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

struct hist_bin {
  double w;
  size_t n;
  hist_bin(): w(0.), n(0) { }
  inline void operator+=(double weight) { w += weight; ++n; }
};

struct ntuple_filler {
  static double weight;
  template <typename T>
  void operator()(T& bin) noexcept { bin += weight; }
};
double ntuple_filler::weight;

using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<ivanp::uniform_axis<double>>>,
  std::vector<hist_bin>,
  ntuple_filler>;

TH1D* root_hist(const hist_t& h, const char* name) {
  TH1D* h1 = new TH1D(name,"",h.axis().nbins(),h.axis().min(),h.axis().max());
  const auto& bins = h.bins();
  size_t n_total = 0;
  for (int i=0, n=bins.size(); i<n; ++i) {
    h1->SetBinContent(i,bins[i].w);
    n_total += bins[i].n;
  }
  h1->SetEntries(n_total);
  return h1;
}

int main(int argc, char* argv[])
{
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm,0.4);
  const double jet_pt_cut = 30.;
  const double jet_eta_cut = 4.4;
  const unsigned need_njets = 2;

  // Define histograms ==============================================
  size_t N = 0, num_events = 0, num_selected = 0;

  hist_t h_H_pT({35,0,700});
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

  std::vector<fastjet::PseudoJet> particles;
  Int_t prev_id = -1;

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
    
    if (prev_id!=*_id) {
      prev_id = *_id;
      N += ( _ncount ? **_ncount : 1);
      ++num_events;
    }

    // Cluster jets -------------------------------------------------
    auto jets = fastjet::ClusterSequence(particles,jet_def)
      .inclusive_jets(jet_pt_cut);
    if (jets.size() < need_njets) continue;
    std::remove_if( jets.begin(), jets.end(),
      [=](const fastjet::PseudoJet& j){ return (j.eta() > jet_eta_cut); });
    if (jets.size() < need_njets) continue;
    std::sort( jets.begin(), jets.end(),
      [](const fastjet::PseudoJet& a, const fastjet::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    // --------------------------------------------------------------

    // Define variables ---------------------------------------------
    ntuple_filler::weight = *_weight;

    const double H_pT = higgs->Pt();
    // --------------------------------------------------------------

    ++num_selected;

    // --------------------------------------------------------------
    h_H_pT(H_pT);
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
  root_hist(h_H_pT,"H_pT");

  fout->Write();

  return 0;
}

