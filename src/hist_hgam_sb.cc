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
#include "nlo_multibin.hh"
#include "Higgs2diphoton.hh"
#include "comprehension.hh"
#include "string_alg.hh"
#include "enum_traits.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;
using ivanp::reserve;

namespace fj = fastjet;
using namespace ivanp::math;

MAKE_ENUM(isp,(any)(gg)(gq)(qq))

isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

template <typename Bin> struct isp_bin {
  std::array<Bin,enum_traits<isp>::size> bins;
  static unsigned i;

  inline Bin& operator->() noexcept { return bins[i]; }
  inline const Bin& operator->() const noexcept { return bins[i]; }
  inline isp_bin& operator+=(double w) noexcept {
    std::get<0>(bins) += w;
    bins[i] += w;
    return *this;
  }
  inline isp_bin& operator+=(const isp_bin& b) noexcept {
    for (auto i=bins.size(); i; ) {
      --i;
      bins[i] += b.bins[i];
    }
    return *this;
  }
};
template <typename Bin> unsigned isp_bin<Bin>::i;

using isp_bin_t = isp_bin<nlo_bin>;
using bin_type = nlo_multibin<isp_bin_t>;

template <typename... Axes>
using hist_t = ivanp::binner<bin_type,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <bool... OF>
using re_hist = ivanp::binner<bin_type, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;

void excl_labels(TH1* h, bool excl) {
  auto* ax = h->GetXaxis();
  for (int i=1, n=h->GetNbinsX(); i<=n; ++i)
    ax->SetBinLabel(i,cat(excl ? "=" : ">=", i-1).c_str());
}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

int main(int argc, char* argv[]) {
  // parse program arguments ========================================
  std::vector<const char*> ntuples, weights;
  const char* output_file_name = nullptr;
  const char* bins_file = STR(PREFIX) "/config/Hjets_ATLAS.bins";
  const char* tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned need_njets = 0;
  bool w_arg = false;
  bool no_photon_cuts = false;

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
    if (!strcmp(argv[i],"--no-photon-cuts")) {
      no_photon_cuts = true;
      continue;
    }

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
  const unsigned njmax = need_njets + 1;

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
  TTreeReaderValue<Int_t> _id1(reader,"id1"), _id2(reader,"id2");

  // handle multiple weights
  std::vector<float_or_double_value_reader> _weights;
  if (!weights_chain) _weights.emplace_back(reader,"weight2");
  else {
    const TObjArray *bb = weights_chain->GetListOfBranches();
    _weights.reserve(bb->GetEntriesFast());
    for (const auto* b : *bb) {
      _weights.emplace_back(reader,static_cast<const TBranch*>(b)->GetName());
    }
  }
  nlo_multibin<>::weights.resize(_weights.size());

  // Define histograms ==============================================
  hist<int> h_N_j_30({need_njets+2u,0,int(need_njets+2)});

#define h_(NAME) re_hist<1> h_##NAME(#NAME,ra[#NAME]);
#define h2_(X1,X2) re_hist<1,0> h_##X1##_##X2(#X1"-"#X2,ra[#X1"_2"],ra[#X2"_2"]);

#define hj_(NAME) auto h_jet_##NAME = reserve<re_hist<1>>(need_njets+1); \
  for (unsigned i=0; i<need_njets+1; ++i) { \
    const auto name = cat(#NAME"_j",i+1); \
    h_jet_##NAME.emplace_back(name,ra[name]); \
  }

  h_(m_yy) h_(pT_yy) h_(cosTS_yy) hj_(pT)

  h2_(m_yy,pT_yy)
  h2_(m_yy,pT_j1)

  h2_(dR_yy,pT_yy)
  h2_(cosTS_yy,pT_yy)

  // ================================================================

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << need_njets
       << "\033[0m or more jets per event\n" << endl;

  Higgs2diphoton Hdecay;
  std::pair<TLorentzVector,TLorentzVector> diphoton;

  size_t ncount = 0, num_events = 0;

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    for (unsigned i=_weights.size(); i!=0; ) { // get weights
      --i; nlo_multibin<>::weights[i] = *_weights[i];
    }

    // Keep track of multi-entry events -----------------------------
    nlo_bin::current_id = *_id;
    const bool new_id = (prev_id != nlo_bin::current_id);
    if (new_id) {
      prev_id = nlo_bin::current_id;
      ncount += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------

    const size_t np = *_nparticle;
    particles.clear();
    particles.reserve(np);
    boost::optional<TLorentzVector> higgs;

    // Read particles -----------------------------------------------
    unsigned n22 = 0;
    for (size_t i=0; i<np; ++i) {
      if (_kf[i] == 25) {
        if (higgs) throw std::runtime_error("more than one Higgs");
        higgs.emplace(_px[i],_py[i],_pz[i],_E[i]);
      } else if (_kf[i] == 22) {
        if (n22>=2) throw std::runtime_error("more than two photons");
        (n22 ? diphoton.second : diphoton.first)
          .SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
        ++n22;
      } else {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }
    if (higgs && n22) throw std::runtime_error("Higgs and photons");
    if (!higgs && n22!=2) throw std::runtime_error("no Higgs or photons");

    isp_bin_t::i = (unsigned)get_isp(*_id1,*_id2);
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

    // Cuts ---------------------------------------------------------
    if (higgs) diphoton = Hdecay(*higgs,new_id);
    else higgs = diphoton.first + diphoton.second;

    auto *A1 = &diphoton.first, *A2 = &diphoton.second;
    double A1_pT = A1->Pt(), A2_pT = A2->Pt();
    if (A1_pT < A2_pT) {
      std::swap(A1,A2);
      std::swap(A1_pT,A2_pT);
    }

    const double H_mass = higgs->M();
    if (!no_photon_cuts) {
      if (A1_pT < 0.35*H_mass) continue;
      if (A2_pT < 0.25*H_mass) continue;
    }

    const double A1_eta = A1->Eta();
    if (!no_photon_cuts)
      if (photon_eta_cut(std::abs(A1_eta))) continue;
    const double A2_eta = A2->Eta();
    if (!no_photon_cuts)
      if (photon_eta_cut(std::abs(A2_eta))) continue;

    h_N_j_30.fill_bin(njets+1); // njets+1 because njets==0 is bin 1

    if (njets < need_njets) continue; // at least needed number of jets

    // Fill histograms ----------------------------------------------
    const double H_pT = higgs->Pt();
    const auto jets_pT = fj_jets | [](const auto& jet){ return jet.pt(); };
    const auto cosTS_yy =
      std::sinh(std::abs(A1_eta-A2_eta)) * A1_pT * A2_pT * 2
      / ( std::sqrt(1.+sq(H_pT/H_mass)) * sq(H_mass) );

    for (unsigned i=0, n=std::min(njets,njmax); i<n; ++i) {
      h_jet_pT[i](jets_pT[i]);
    }

    h_m_yy(H_mass);
    h_pT_yy(H_pT);
    h_cosTS_yy(cosTS_yy);

    h_m_yy_pT_yy(H_mass,H_pT);
    h_m_yy_pT_j1(H_mass,jets_pT[0]);

    h_dR_yy_pT_yy(diphoton.first.DeltaR(diphoton.second),H_pT);

    h_cosTS_yy_pT_yy(cosTS_yy,H_pT);

  } // END EVENT LOOP
  // ================================================================

  cout << "Processed events: " << num_events << endl;
  cout << "ncount: " << ncount << '\n' << endl;

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(output_file_name,"recreate");
  if (fout->IsZombie()) return 1;
  TDirectory *dir = fout.get();

  auto h_N_j_30_integrated = h_N_j_30;
  h_N_j_30_integrated.integrate_left();

  // write root historgrams
  isp_bin_t::i = 0;
  for (const char* isp_str : enum_traits<isp>::all_str()) {
    if (isp_bin_t::i) dir = dir->mkdir(isp_str);
    bin_type::wi = 0;
    for (const auto& _w : _weights) {
      dir = dir->mkdir(cat(_w.GetBranchName(),"_Jet",
          jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
        : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
        : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
        : "", std::setprecision(2), jet_def.R()*10.
      ).c_str());
      cout << dir->GetName() << endl;
      dir->cd();

      using ivanp::root::to_root;
      using ivanp::root::slice_to_root;

      auto* h_N_j_30_excl = to_root(h_N_j_30,"N_j_30_excl");
      excl_labels(h_N_j_30_excl,true);
      auto* h_N_j_30_incl = to_root(h_N_j_30_integrated,"N_j_30_incl");
      h_N_j_30_incl->SetEntries( h_N_j_30_excl->GetEntries() );
      excl_labels(h_N_j_30_incl,false);

      for (auto& h : re_hist<1>::all) to_root(*h,h.name);
      for (auto& h : re_hist<1,0>::all) {
        const auto vars = ivanp::rsplit<1>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1]);
      }

      ++bin_type::wi;
      dir = dir->GetMotherDir();
    }
    if (isp_bin_t::i) dir = dir->GetMotherDir();
    ++isp_bin_t::i;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_events);
  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}
