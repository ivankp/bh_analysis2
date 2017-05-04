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
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "re_axes.hh"
#include "parse_args.hh"
#include "lo_multibin.hh"
#include "lo_profile_multibin.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

// #define NLO

using std::cout;
using std::cerr;
using std::endl;

using boost::optional;

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

template <typename... Axes>
using hist_t = ivanp::binner<lo_multibin,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <bool... OF>
using re_hist = ivanp::binner<lo_multibin, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;
template <bool... OF>
using re_prof = ivanp::binner<lo_profile_multibin, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;

int main(int argc, char* argv[]) {
  // parse program arguments ========================================
  std::vector<const char*> ntuples, weights;
  const char* output_file_name = nullptr;
  const char* bins_file = STR(PREFIX) "/config/Hjets_isolation.bins";
  const char* tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned need_njets = 0;
  bool w_arg = false;

  for (int i=1; i<argc; ++i) {
    using namespace parse;
    if (!std::strcmp(argv[i],"w")) { w_arg = true; continue; }
    if (!std::strcmp(argv[i],"bh")) { w_arg = false; continue; }
    if (output_root_file_name(argv[i],output_file_name)) continue;
    if (root_file_name(argv[i], w_arg ? weights : ntuples )) continue;
    if (num_jets(argv[i],need_njets)) continue;
    if (bins_file_name(argv[i],bins_file)) continue;
    if (jetdef(argv[i],jet_def)) continue;
    if (parse::tree_name(argv[i],tree_name)) continue;

    cerr << "\033[31mUnexpected argument:\033[0m " << argv[i] << endl;
    return 1;
  }
  if (!ntuples.size()) {
    cerr << "\033[31mNo input root files provided\033[0m" << endl;
    return 1;
  }
  if (!output_file_name) {
    cerr << "\033[31mNo output root file provided\033[0m" << endl;
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
#ifdef NLO
  const unsigned njmax = need_njets + 1;
#else
  const unsigned njmax = need_njets;
#endif

  cout << "\033[36mBinning\033[0m: " << bins_file << '\n' << endl;
  re_axes ra(bins_file);

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (const char* name : ntuples) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  optional<TChain> weights_chain;
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

  optional<TTreeReaderValue<Int_t>> _ncount;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
  }

  // handle multiple weights
  std::vector<float_or_double_value_reader> _weights;
  if (!weights_chain) _weights.emplace_back(reader,"weight");
  else {
    const TObjArray *bb = weights_chain->GetListOfBranches();
    _weights.reserve(bb->GetEntriesFast());
    for (const auto* b : *bb) {
      _weights.emplace_back(reader,static_cast<const TBranch*>(b)->GetName());
    }
  }
  lo_multibin::weights.resize(_weights.size());

  // Define histograms ==============================================
  hist<int> h_Njets({njmax+1u,0,int(njmax+1)});

#define a_(name) auto a_##name = ra[#name];
#define h_(name) re_hist<1> h_##name(#name,ra[#name]);
#define p_(name) re_prof<1> p_##name(#name,ra[#name]);

  a_(dR)

  h_(H_jet_dR_min)
  auto h_H_jet_dR = reserve<re_hist<1>>(njmax);
  for (unsigned i=0; i<njmax; ++i) {
    h_H_jet_dR.emplace_back(cat("H_jet",i+1,"_dR"),a_dR);
  }

  a_(y) a_(phi)

  p_(HT) p_(H_pT) p_(H_y) p_(H_eta) p_(H_phi)

  auto p_jet_pT   = reserve<re_prof<1>>(njmax);
  auto p_jet_y    = reserve<re_prof<1>>(njmax);
  auto p_jet_eta  = reserve<re_prof<1>>(njmax);
  auto p_jet_phi  = reserve<re_prof<1>>(njmax);

  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_pT");
    p_jet_pT.emplace_back(name,ra[name]);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_y");
    p_jet_y.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_eta");
    p_jet_eta.emplace_back(name,a_y);
  }
  for (unsigned i=0; i<njmax; ++i) {
    auto name = cat("jet",i+1,"_phi");
    p_jet_phi.emplace_back(name,a_phi);
  }

#define jjpT_p_(name) \
  optional<re_prof<1>> p_jjpT_##name; \
  if (need_njets >= 2) { \
    const char* _name = need_njets > 2 ? "jjpT_"#name : "jj_"#name; \
    p_jjpT_##name.emplace( _name, ra[_name] ); \
  }
#define jjfb_p_(name) \
  optional<re_prof<1>> p_jjfb_##name; \
  if (need_njets > 2) p_jjfb_##name.emplace("jjfb_"#name,ra["jjfb_"#name]);

  jjpT_p_(dpT )  jjfb_p_(dpT )
  jjpT_p_(dy  )  jjfb_p_(dy  )
  jjpT_p_(deta)  jjfb_p_(deta)
  jjpT_p_(dphi)  jjfb_p_(dphi)
  jjpT_p_(mass)  jjfb_p_(mass)

  p_(Hjets_mass)

  // ================================================================

  std::vector<fj::PseudoJet> particles;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n" << endl;

  size_t ncount = 0, num_selected = 0;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    for (unsigned i=_weights.size(); i!=0; ) { // get weights
      --i; multibin::weights[i] = *_weights[i];
    }

    const size_t np = *_nparticle;
    particles.clear();
    particles.reserve(np);
    optional<TLorentzVector> higgs;

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
    TLorentzVector Hjets = *higgs;
    for (const auto& j : jets) {
      HT += j.pT;
      Hjets += j.p;
    }

    const double H_y = higgs->Rapidity();
    const double H_phi = higgs->Phi();
    // --------------------------------------------------------------

    ++num_selected;

    double min_dR = 1e55;
    for (unsigned i=0, n=std::min(njets,njmax); i<n; ++i) {
      const double dR = deltaR(H_y,jets[i].y,H_phi,jets[i].phi);
      h_H_jet_dR[i](dR);
      if (dR < min_dR) min_dR = dR;
    }
    h_H_jet_dR_min(min_dR);

    // Fill histograms ----------------------------------------------
    p_HT(HT,min_dR);

    p_H_pT(H_pT,min_dR);
    p_H_y(H_y,min_dR);
    p_H_eta(higgs->Eta(),min_dR);
    p_H_phi(H_phi,min_dR);

    // jet histograms
    for (unsigned i=0, n=std::min(njets,njmax); i<n; ++i) {
      p_jet_pT  [i](jets[i].pT  ,min_dR);
      p_jet_y   [i](jets[i].y   ,min_dR);
      p_jet_eta [i](jets[i].eta ,min_dR);
      p_jet_phi [i](jets[i].phi ,min_dR);
    }

    p_Hjets_mass(Hjets.M(),min_dR);

    if (njets<2) continue; // 222222222222222222222222222222222222222

    // Jet pair with highest pT .....................................
    const dijet jjpT(jets[0],jets[1]);

    (*p_jjpT_dpT )(jjpT.dpT ,min_dR);
    (*p_jjpT_dy  )(jjpT.dy  ,min_dR);
    (*p_jjpT_deta)(jjpT.deta,min_dR);
    (*p_jjpT_dphi)(jjpT.dphi,min_dR);
    (*p_jjpT_mass)(jjpT.mass,min_dR);
    // ..............................................................

    if (njets<3) continue; // 333333333333333333333333333333333333333

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

    (*p_jjfb_dpT )(jjfb.dpT ,min_dR);
    (*p_jjfb_dy  )(jjfb.dy  ,min_dR);
    (*p_jjfb_deta)(jjfb.deta,min_dR);
    (*p_jjfb_dphi)(jjfb.dphi,min_dR);
    (*p_jjfb_mass)(jjfb.mass,min_dR);
    // ..............................................................

  } // END EVENT LOOP
  // ================================================================

  cout << "Selected entries: " << num_selected << endl;
  cout << "ncount: " << ncount << '\n' << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(output_file_name,"recreate");
  if (fout->IsZombie()) return 1;

  // write root historgrams
  lo_multibin::wi = 0;
  for (const auto& _w : _weights) {
    auto* dir = fout->mkdir(cat(_w.GetBranchName(),"_Jet",
        jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
      : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
      : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
      : "", std::setprecision(2), jet_def.R()*10.
    ).c_str());
    cout << dir->GetName() << endl;
    dir->cd();

    using ivanp::root::to_root;
    using ivanp::root::slice_to_root;

    auto* h_Njets_excl = to_root(h_Njets,"jets_N_excl");
    h_Njets.integrate_left();
    auto* h_Njets_incl = to_root(h_Njets,"jets_N_incl");
    h_Njets_incl->SetEntries( h_Njets_excl->GetEntries() );

    for (auto& h : re_hist<1>::all) to_root(*h,h.name);

    dir->mkdir("profiles")->cd();
    for (auto& p : re_prof<1>::all) to_root(*p,p.name);
    dir->cd();

    ++lo_multibin::wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_selected);
  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}

