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
#include "bin_defs.hh"
#include "re_axes.hh"
#include "parse_args.hh"
#include "Higgs2diphoton.hh"
#include "comprehension.hh"
#include "string_alg.hh"
#include "enum_traits.hh"
#include "category_bin.hh"

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

MAKE_ENUM(isp,(all)(gg)(gq)(qq))

isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

MAKE_ENUM(photon_cuts,(all)(with_photon_cuts))

auto& wi = multiweight_bin_base::wi;

using cat_bin = category_bin<nlo_bin,photon_cuts,isp>;

using bin_t = multiweight_bin<cat_bin>;
template <bool... OF>
using hist = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<typename re_axes::axis_type,OF,OF>...> >;

void excl_labels(TH1* h, bool excl) {
  auto* ax = h->GetXaxis();
  for (int i=1, n=h->GetNbinsX(); i<=n; ++i)
    ax->SetBinLabel(i,cat(excl ? "=" : ">=", i-1).c_str());
}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

TLorentzVector operator+(const TLorentzVector& a, const fastjet::PseudoJet& b) {
  return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] };
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
  bin_t::weights.resize(_weights.size());

  // Define histograms ==============================================
  ivanp::binner<bin_t, std::tuple<ivanp::axis_spec<
      ivanp::uniform_axis<int>
    >>> h_N_j_30({need_njets+2u,0,int(need_njets+2)});

#define h_(NAME) hist<1> h_##NAME(#NAME,ra[#NAME]);
#define h2_(X1,X2) hist<1,0> h_##X1##_##X2(#X1"-"#X2,ra[#X1"_2"],ra[#X2"_2"]);

#define hj_(NAME) auto h_jet_##NAME = reserve<hist<1>>(need_njets+1); \
  for (unsigned i=0; i<need_njets+1; ++i) { \
    const auto name = cat(#NAME"_j",i+1); \
    h_jet_##NAME.emplace_back(name,ra[name]); \
  }

  h_(m_yy) h_(pT_yy)
  h_(pT_yy_105_160) h_(pT_yy_121_129)
  hj_(pT)
  h_(cosTS_yy) h_(cos_y1) h_(cos_y2)

  h2_(m_yy,pT_yy)
  h2_(m_yy,pT_j1)

  h2_(dR_yy,pT_yy)
  h2_(cosTS_yy,pT_yy)
  h2_(cos_y1,pT_yy)

  h2_(m_yyj,eta_yy) h2_(m_yyj,eta_j1)

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
      --i; bin_t::weights[i] = *_weights[i];
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

    cat_bin::id<isp>() = (unsigned)get_isp(*_id1,*_id2);
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

    TLorentzVector *A1 = &diphoton.first, *A2 = &diphoton.second;
    double A1_pT = A1->Pt(), A2_pT = A2->Pt();
    if (A1_pT < A2_pT) {
      std::swap(A1,A2);
      std::swap(A1_pT,A2_pT);
    }

    bool passes_photon_cuts = true;

    const double H_mass = higgs->M();
    if (A1_pT < 0.35*125.) passes_photon_cuts = false;
    if (A2_pT < 0.25*125.) passes_photon_cuts = false;

    const double A1_eta = A1->Eta();
    const double A2_eta = A2->Eta();
    if (photon_eta_cut(std::abs(A1_eta))) passes_photon_cuts = false;
    if (photon_eta_cut(std::abs(A2_eta))) passes_photon_cuts = false;

    cat_bin::id<photon_cuts>() = ( passes_photon_cuts ? 1 : 0 );

    h_N_j_30.fill_bin(njets+1); // njets+1 because njets==0 is bin 1

    if (njets < need_njets) continue; // at least needed number of jets

    // Fill histograms ----------------------------------------------
    const double H_pT = higgs->Pt();
    const auto jets_pT = fj_jets | [](const auto& jet){ return jet.pt(); };

    const double cosTS_yy =
      std::sinh(std::abs(A1_eta-A2_eta)) * A1_pT * A2_pT * 2
      / ( std::sqrt(1.+sq(H_pT/H_mass)) * sq(H_mass) );

    TLorentzVector AAf = *higgs, A1f = *A1, A2f = *A2;
    const auto boost = -AAf.BoostVector();
    AAf.Boost(boost);
    A1f.Boost(boost);
    A2f.Boost(boost);
    const double A1_cos = std::abs(A1f.CosTheta());
    const double A2_cos = std::abs(A2f.CosTheta());

    for (unsigned i=0, n=std::min(njets,njmax); i<n; ++i) {
      h_jet_pT[i](jets_pT[i]);
    }

    h_m_yy(H_mass);
    h_pT_yy(H_pT);

    if (105.<H_mass && H_mass<160.) {
      h_pT_yy_105_160(H_pT);
      if (121.<H_mass && H_mass<129.)
        h_pT_yy_121_129(H_pT);
    }

    h_cosTS_yy(cosTS_yy);
    h_cos_y1(A1_cos);
    h_cos_y2(A2_cos);

    h_m_yy_pT_yy(H_mass,H_pT);
    h_m_yy_pT_j1(H_mass,jets_pT[0]);

    h_dR_yy_pT_yy(diphoton.first.DeltaR(diphoton.second),H_pT);

    h_cosTS_yy_pT_yy(cosTS_yy,H_pT);
    h_cos_y1_pT_yy(A1_cos,H_pT);

    const auto Hj = *higgs + fj_jets[0];
    const double Hj_mass = Hj.M();
    const double H_eta = higgs->Eta();
    const double j1_eta = fj_jets[0].eta();

    h_m_yyj_eta_yy(Hj_mass,H_eta);
    h_m_yyj_eta_j1(Hj_mass,j1_eta);

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

#define CATEGORY_TOP(NAME) \
  cat_bin::id<NAME>() = 0; \
  for (const char* NAME##_str : enum_traits<NAME>::all_str()) { \
    dir = dir->mkdir(NAME##_str);

#define CATEGORY_BOT(NAME) \
    dir = dir->GetMotherDir(); \
    ++cat_bin::id<NAME>(); \
  }

  // write root historgrams
  wi = 0;
  for (const auto& _w : _weights) {
    dir = dir->mkdir(cat(_w.GetBranchName(),"_Jet",
        jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
      : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
      : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
      : "", std::setprecision(2), jet_def.R()*10.
    ).c_str());
    cout << dir->GetName() << endl;

    CATEGORY_TOP(photon_cuts)
    CATEGORY_TOP(isp)

      dir->cd();

      using ivanp::root::to_root;
      using ivanp::root::slice_to_root;

      auto* h_N_j_30_excl = to_root(h_N_j_30,"N_j_30_excl");
      excl_labels(h_N_j_30_excl,true);
      auto* h_N_j_30_incl = to_root(h_N_j_30_integrated,"N_j_30_incl");
      h_N_j_30_incl->SetEntries( h_N_j_30_excl->GetEntries() );
      excl_labels(h_N_j_30_incl,false);

      for (auto& h : hist<1>::all) to_root(*h,h.name);
      for (auto& h : hist<1,0>::all) {
        const auto vars = ivanp::rsplit<1>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1]);
      }

    CATEGORY_BOT(isp)
    CATEGORY_BOT(photon_cuts)

    dir = dir->GetMotherDir();
    ++wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_events);
  fout->Write();
  cout << "\n\033[32mOutput\033[0m: " << fout->GetName() << endl;

  return 0;
}

