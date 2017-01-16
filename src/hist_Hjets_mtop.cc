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
#include "re_axes.hh"
#include "timed_counter.hh"
#include "catstr.hh"

#define test(var) \
  cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << endl;

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;

using p4_t = Float_t;

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

struct dijet {
  TLorentzVector p;
  double dpT, dy, deta, dphi, mass;
  dijet(const Jet& j1, const Jet& j2)
  : p(j1.p+j2.p),
    dpT (std::abs(j1.pT -j2.pT )),
    dy  (std::abs(j1.y  -j2.y  )),
    deta(std::abs(j1.eta-j2.eta)),
    dphi(ivanp::math::dphi(j1.phi,j2.phi)),
    mass(p.M()) {}
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
};
double hist_bin::weight;

template <typename T>
using axis = ivanp::uniform_axis<T>;
template <typename... Axes>
using hist_t = ivanp::binner<hist_bin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<axis<T>>;

using re_axis = typename re_axes::axis_type;
template <size_t N>
using re_hist = ivanp::binner<hist_bin,
  ivanp::tuple_of_same_t<ivanp::axis_spec<re_axis>,N>>;

template <typename A>
TH1D* root_hist(const hist_t<A>& h, const std::string& name) {
  TH1D* hr = new TH1D(name.c_str(),"",
                      h.axis().nbins(),h.axis().min(),h.axis().max());
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
  constexpr unsigned need_njets = 2;

  // Define histograms ==============================================
  size_t N = 0, num_events = 0, num_selected = 0;

  axis<int> a_Njets(need_njets+2,0,need_njets+2);
  hist<int> h_Njets_incl(a_Njets), h_Njets_excl(a_Njets);

  re_axes ra("binning.txt");
#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);

  a_(y) a_(absy) a_(phi)

  h_(HT) h_(H_pT) h_(H_y) h_(H_eta) h_(H_phi) h_(H_mass)

  std::vector<re_hist<1>> h_jet_pT, h_jet_y, h_jet_eta, h_jet_phi, h_jet_mass;
    h_jet_pT.reserve(need_njets+1);
    h_jet_y.reserve(need_njets+1);
    h_jet_eta.reserve(need_njets+1);
    h_jet_phi.reserve(need_njets+1);
    h_jet_mass.reserve(need_njets+1);
  for (unsigned i=0; i<need_njets+1; ++i) {
    static std::string name;
    name = cat("jet",i+1,"_pT");
    h_jet_pT.emplace_back(name,ra[name]);
    name = cat("jet",i+1,"_y");
    h_jet_y.emplace_back(name,a_y);
    name = cat("jet",i+1,"_eta");
    h_jet_eta.emplace_back(name,a_y);
    name = cat("jet",i+1,"_phi");
    h_jet_phi.emplace_back(name,a_phi);
    name = cat("jet",i+1,"_mass");
    h_jet_mass.emplace_back(name,ra[name]);
  }

  h_(jjpT_dpT )  h_(jjfb_dpT )
  h_(jjpT_dy  )  h_(jjfb_dy  )
  h_(jjpT_deta)  h_(jjfb_deta)
  h_(jjpT_dphi)  h_(jjfb_dphi)
  h_(jjpT_mass)  h_(jjfb_mass)

  // histograms for mtop study
  a_(_x) a_(_HT) a_(_maxdy) a_(_pT) a_(_x2)

  re_hist<3>
    h_xH_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x1_HT_maxdy(a__x, a__HT, a__maxdy),
    h_x2_HT_maxdy(a__x, a__HT, a__maxdy);

  // p1's pT in bins of p2's x
  std::array<std::array< re_hist<2>, 3>,3> h_p1pT_p2x;
  for (auto& a : h_p1pT_p2x) for (auto& h : a) h = {a__pT,a__x2};

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
    hist_bin::weight = *_weight; // Read weight

    // Keep track of multi-entry events -----------------------------
    curr_id = *_id;
    if (prev_id!=curr_id) {
      prev_id = curr_id;
      N += ( _ncount ? **_ncount : 1);
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

    // Fill Njets histograms before cuts ----------------------------
    h_Njets_excl.fill_bin(njets+1); // nj+1 because nj==0 is bin 1
    for (unsigned j=1; j<=njets+1; ++j) h_Njets_incl.fill_bin(j);
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
    // --------------------------------------------------------------

    ++num_selected;

    // Fill histograms ----------------------------------------------
    h_HT(HT);

    h_H_pT(H_pT);
    h_H_y(H_y);
    h_H_eta(higgs->Eta());
    h_H_phi(higgs->Phi());
    h_H_mass(higgs->M());

    // jet histograms
    for (unsigned i=0, n=std::min(njets,need_njets+1); i<n; ++i) {
      h_jet_pT  [i](jets[i].pT  );
      h_jet_y   [i](jets[i].y   );
      h_jet_eta [i](jets[i].eta );
      h_jet_phi [i](jets[i].phi );
      h_jet_mass[i](jets[i].mass);
    }

    // find maximum rapidity separation in the event
    double max_dy = abs(jets[0].y-H_y);
    for (unsigned i=1; i<njets; ++i) {
      for (unsigned j=0; j<i; ++j) {
        larger(max_dy,abs(jets[i].y-jets[j].y));
      }
      larger(max_dy,abs(jets[i].y-H_y));
    }

    if (njets<2) continue; // 222222222222222222222222222222222222222

    // HT fraction x in bins on HT and max rapidity separation
    const double xH = H_pT/HT, x1 = jets[0].pT/HT, x2 = jets[1].pT/HT;
    h_xH_HT_maxdy(xH, HT, max_dy);
    h_x1_HT_maxdy(x1, HT, max_dy);
    h_x2_HT_maxdy(x2, HT, max_dy);

    // p1's pT in bins of p2's x
    h_p1pT_p2x[0][0](H_pT, xH);
    h_p1pT_p2x[1][0](jets[0].pT, xH);
    h_p1pT_p2x[2][0](jets[1].pT, xH);
    h_p1pT_p2x[0][1](H_pT, x1);
    h_p1pT_p2x[1][1](jets[0].pT, x1);
    h_p1pT_p2x[2][1](jets[1].pT, x1);
    h_p1pT_p2x[0][2](H_pT, x2);
    h_p1pT_p2x[1][2](jets[0].pT, x2);
    h_p1pT_p2x[2][2](jets[1].pT, x2);

    // Jet pair with highest pT .....................................
    const dijet jjpT(jets[0],jets[1]);

    h_jjpT_dpT (jjpT.dpT );
    h_jjpT_dy  (jjpT.dy  );
    h_jjpT_deta(jjpT.deta);
    h_jjpT_dphi(jjpT.dphi);
    h_jjpT_mass(jjpT.mass);
    // ..............................................................

    // Two most forward-backward jets ...............................
    unsigned jb = 0, jf = 0;
    if (njets==2) { jf = 1;
    } else {
      for (unsigned j=1; j<njets; ++j) {
        if (jets[j].y < jets[jb].y) jb = j;
        if (jets[j].y > jets[jf].y) jf = j;
      }
    }
    const dijet& jjfb = ( (jb==0 && jf==1) || (jb==1 && jf==0)
      ? jjpT : dijet(jets[jb],jets[jf]) // possibly reuse jjpT
    );

    h_jjfb_dpT (jjfb.dpT );
    h_jjfb_dy  (jjfb.dy  );
    h_jjfb_deta(jjfb.deta);
    h_jjfb_dphi(jjfb.dphi);
    h_jjfb_mass(jjfb.mass);
    // ..............................................................

    // --------------------------------------------------------------
  } // END EVENT LOOP
  // ****************************************************************

  test(N)

  cout << "Selected entries: " << num_selected << endl;
  cout << "Processed events: " << num_events << endl;

  // Open input ntuple root file
  auto fout = std::make_unique<TFile>(argv[2],"recreate");
  if (fout->IsZombie()) return 1;

  fout->mkdir("weight2_JetAntiKt4")->cd();

  // write root historgrams
#define rh(name) root_hist(h_##name,#name);
  rh(Njets_incl)
  rh(Njets_excl)

  rh(HT)

  rh(H_pT)
  rh(H_y)
  rh(H_eta)
  rh(H_phi)
  rh(H_mass)

  for (unsigned i=0; i<=need_njets; ++i)
    root_hist(h_jet_pT  [i],cat("jet",i+1,"_pT"));
  for (unsigned i=0; i<=need_njets; ++i)
    root_hist(h_jet_y   [i],cat("jet",i+1,"_y"));
  for (unsigned i=0; i<=need_njets; ++i)
    root_hist(h_jet_eta [i],cat("jet",i+1,"_eta"));
  for (unsigned i=0; i<=need_njets; ++i)
    root_hist(h_jet_phi [i],cat("jet",i+1,"_phi"));
  for (unsigned i=0; i<=need_njets; ++i)
    root_hist(h_jet_mass[i],cat("jet",i+1,"_mass"));

  if (need_njets >= 2) {
    rh(jjpT_dpT ) rh(jjfb_dpT )
    rh(jjpT_dy  ) rh(jjfb_dy  )
    rh(jjpT_deta) rh(jjfb_deta)
    rh(jjpT_dphi) rh(jjfb_dphi)
    rh(jjpT_mass) rh(jjfb_mass)
  }

  fout->cd();
  (new TH1D("N","N",1,0,1))->SetBinContent(1,N);
  fout->Write();

  return 0;
}

