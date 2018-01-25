// Written by Ivan Pogrebnyak

#include <iostream>
#include <algorithm>

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
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "bin_defs.hh"
#include "re_axes.hh"
#include "Higgs2diphoton.hh"
#include "comprehension.hh"
#include "program_options.hh"
#include "parse_args/jetdef.hh"
#include "tc_msg.hh"
#include "string_alg.hh"

#define TEST(VAR) \
  std::cout << tc::cyan << #VAR << tc::reset << " = " << VAR << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
namespace tc = termcolor;
namespace fj = fastjet;
using namespace ivanp;
using namespace ivanp::math;

template <typename... Args>
inline void write(const char* name, Args&&... args) {
  TNamed(name,cat(args...).c_str()).Write();
}

using bin_t = multiweight_bin<nlo_bin>;
template <typename... Axes>
using hist_t = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<Axes>...>>;
template <typename T>
using hist = hist_t<ivanp::uniform_axis<T>>;

using re_axis = typename re_axes::axis_type;
template <bool... OF>
using re_hist = ivanp::binner<bin_t, std::tuple<
  ivanp::axis_spec<re_axis,OF,OF>...>>;

void excl_labels(TH1* h, bool excl) {
  auto* ax = h->GetXaxis();
  for (int i=1, n=h->GetNbinsX(); i<=n; ++i)
    ax->SetBinLabel(i,cat(excl ? "=" : ">=", i-1).c_str());
}

template <typename T>
inline TLorentzVector operator+(const TLorentzVector& a, const T& b) {
  return { a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3] };
}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ntuples, weights;
  const char* ofname = nullptr;
  const char* bfname = STR(PREFIX) "/config/Hjets_ang.bins";
  const char* tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned njets_expected;
  double jet_pt_cut = 30., jet_eta_cut = 4.4;
  bool no_photon_cuts = false;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ntuples,'i',"input ROOT BH ntuples",req(),pos())
      (weights,'w',"input ROOT weight ntuples")
      (ofname,'o',"output file name",req())
      (bfname,'b',cat("bins file name [",bfname,"]"))
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (njets_expected,'j',"expected number of jets",req())
      (jet_def,{},"jet clustering algorithm [AntiKt4]",
       parse::jetdef,default_init(fj::antikt_algorithm,0.4))
      (jet_pt_cut,"--jet-pt-cut",cat('[',jet_pt_cut,"] GeV"))
      (jet_eta_cut,"--jet-eta-cut",cat('[',jet_eta_cut,']'))
      (no_photon_cuts,"--no-photon-cuts")
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  // const unsigned njmax = njets_expected + 1;

  info("Binning", bfname); cout << '\n';
  re_axes ra(bfname);

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  info("Input ntuples");
  for (const char* name : ntuples) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  TChain *weights_chain = nullptr;
  if (weights.size()) {
    info("Input weights");
    weights_chain = new TChain("weights");
    for (const char* name : weights) {
      if (!weights_chain->Add(name,0)) return 1;
      cout << "  " << name << endl;
    }
    cout << endl;
  }

  // Set up branches for reading
  if (weights_chain) chain.AddFriend(weights_chain);
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
    if (_ncount) break;
  }

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
#define h_(NAME) re_hist<1> h_##NAME(#NAME,ra[#NAME]);

  h_(H_cosTheta)
  h_(Hj_mass)

  re_hist<1,0> h_H_absCosTheta_Hj_mass(
    "H_absCosTheta-Hj_mass", ra["H_absCosTheta_2"], ra["Hj_mass_2"]);

  // ================================================================

  std::vector<fj::PseudoJet> particles;
  particles.reserve(njets_expected+1);
  Int_t prev_id = -1;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting " << njets_expected
       << " or more jets per event\n" << endl;

  TLorentzVector Higgs;
  Higgs2diphoton Hdecay;

  size_t ncount = 0, num_events = 0;

  // LOOP ===========================================================
  using cnt = ivanp::timed_counter<Long64_t>;
  for (cnt ent(reader.GetEntries(true)); reader.Next(); ++ent) {
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

    // Read particles -----------------------------------------------
    const size_t np = *_nparticle;
    particles.clear();
    for (size_t i=0; i<np; ++i) {
      if (_kf[i] == 25) {
        Higgs = {_px[i],_py[i],_pz[i],_E[i]};
      } else {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }
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

    if (njets < njets_expected) continue; // require expected number of jets

    // Photons ------------------------------------------------------
    if (!no_photon_cuts) {
      const double H_mass = Higgs.M();

      const auto photons = Hdecay(Higgs);
      auto *A1 = &photons.first, *A2 = &photons.second;
      double A1_pT = A1->Pt(), A2_pT = A2->Pt();
      if (A1_pT < A2_pT) {
        std::swap(A1,A2);
        std::swap(A1_pT,A2_pT);
      }

      if (A1_pT < 0.35*H_mass) continue;
      if (A2_pT < 0.25*H_mass) continue;

      if (photon_eta_cut(std::abs(A1->Eta()))) continue;
      if (photon_eta_cut(std::abs(A2->Eta()))) continue;
    }
    // --------------------------------------------------------------

    const TLorentzVector jet1(
      fj_jets.front().px(),
      fj_jets.front().py(),
      fj_jets.front().pz(),
      fj_jets.front().e()
    );
    const auto Q = Higgs + jet1;
    const double Q2 = Q*Q;
    const double M = sqrt(Q2);

    const TLorentzVector Z(0,0,Q.E(),Q.Pz());
    const auto ell = ((Q*jet1)/Q2)*Higgs - ((Q*Higgs)/Q2)*jet1;

    const double cos_theta = (ell*Z) / sqrt(sq(ell)*sq(Z));

    h_H_cosTheta(cos_theta);
    h_Hj_mass(M);

    h_H_absCosTheta_Hj_mass(std::abs(cos_theta),M);

  } // END EVENT LOOP
  // ================================================================
  cout.imbue(std::locale(""));
  info("Processed events", num_events);
  info("ncount", ncount); cout << '\n';

  // Open output root file for histograms
  TFile fout(ofname,"recreate","",109);
  if (fout.IsZombie()) return 1;

  // write root historgrams
  bin_t::wi = 0;
  auto mkdir = [&fout, jet_alg = cat("_Jet",
        jet_def.jet_algorithm() == fj::antikt_algorithm ? "AntiKt"
      : jet_def.jet_algorithm() == fj::kt_algorithm ? "Kt"
      : jet_def.jet_algorithm() == fj::cambridge_algorithm ? "CA"
      : "", std::setprecision(2), jet_def.R()*10.
    )](const char* branch_name){
      return fout.mkdir((branch_name+jet_alg).c_str());
    };
  info("Directories");
  for (const auto& _w : _weights) {
    auto* dir = mkdir(_w.GetBranchName());
    cout << "  " << dir->GetName() << endl;
    dir->cd();

    using ivanp::root::to_root;
    using ivanp::root::slice_to_root;

    for (auto& h : re_hist<1>::all) to_root(*h,h.name);

    for (auto& h : re_hist<1,0>::all) {
      const auto vars = rsplit<1>(h.name,'-');
      slice_to_root(*h,vars[0],vars[1]);
    }

    ++bin_t::wi;
  }

  fout.cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount);
  h_N->SetEntries(num_events);
  fout.Write();
  info("Output file", fout.GetName());

  return 0;
}
