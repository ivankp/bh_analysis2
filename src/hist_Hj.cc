#include <iostream>
#include <numeric>

#include <boost/optional.hpp>

#include <TFile.h>
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
#include "float_or_double_reader.hh"
#include "binner_root.hh"
#include "bin_defs.hh"
#include "category_bin.hh"
#include "re_axes.hh"
#include "Higgs2diphoton.hh"
#include "comprehension.hh"
#include "string_alg.hh"
#include "program_options.hh"
#include "parse_args/jetdef.hh"
#include "tc_msg.hh"

#define TEST(VAR) \
  std::cout << tc::cyan << #VAR << tc::reset << " = " << VAR << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;
using ivanp::reserve;

namespace fj = fastjet;
namespace tc = termcolor;
using namespace ivanp::math;

template <typename Name, typename Title>
inline void note(Name&& name, Title&& title) {
  TNamed(std::forward<Name>(name),std::forward<Title>(title)).Write();
}

MAKE_ENUM(isp,(all)(gg)(gq)(qq))

isp get_isp(Int_t id1, Int_t id2) noexcept {
  const bool g1 = (id1 == 21), g2 = (id2 == 21);
  if (g1 == g2) return g1 ? isp::gg : isp::qq;
  else return isp::gq;
}

MAKE_ENUM(photon_cuts,(all)(with_photon_cuts))

#define IMPL_GLOBAL
#include STR(IMPL)
#undef IMPL_GLOBAL

#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/reverse.hpp>

#ifndef CATEGORIES
#define CATEGORIES (photon_cuts)(isp)
#endif

using cat_bin = category_bin<nlo_bin,BOOST_PP_SEQ_ENUM(CATEGORIES)>;

using bin_t = multiweight_bin<cat_bin>;
template <bool... OF>
using hist = ivanp::binner<bin_t,
  std::tuple<ivanp::axis_spec<typename re_axes::axis_type,OF,OF>...> >;

auto& wi = multiweight_bin_base::wi;

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
  const char *ofname = nullptr, *bfname = nullptr;
  const char *tree_name = "t3";
  fj::JetDefinition jet_def;
  unsigned njets_expected;
  double jet_pt_cut = 30., jet_eta_cut = 4.4;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ntuples,'i',"input ROOT BH ntuples",req(),pos())
      (weights,'w',"input ROOT weight ntuples")
      (ofname,'o',"output file name",req())
#ifndef BINS_FILE
      (bfname,'b',"bins file name",req())
#else
      (bfname,'b',"bins file name [" BINS_FILE "]",
       default_init(STR(PREFIX) "/" BINS_FILE))
#endif
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (njets_expected,'j',"expected number of jets",req())
      (jet_def,{},"jet clustering algorithm [AntiKt4]",
       parse::jetdef,default_init(fj::antikt_algorithm,0.4))
      (jet_pt_cut,"--jet-pt-cut",cat('[',jet_pt_cut,"] GeV"))
      (jet_eta_cut,"--jet-eta-cut",cat('[',jet_eta_cut,']'))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }
  // ================================================================

  cout << "\033[36mBinning\033[0m: " << bfname << '\n' << endl;
  re_axes ra(bfname);

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
    if (!strcmp(bo->GetName(),"ncount")) {
      _ncount.emplace(reader,"ncount"); break;
    }
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
      ivanp::uniform_axis<int>,0,1
    >>> h_Njets({njets_expected+2u,0,int(njets_expected+2)});

#define H_MACRO(_1,_2,_3,NAME,...) NAME
#define h_(...) H_MACRO(__VA_ARGS__, h3_, h2_, h1_)(__VA_ARGS__)

#define h1_(X) hist<1> h_##X(#X,ra[#X]);
#define h2_(X1,X2) hist<1,0> h_##X1##_##X2(#X1"-"#X2,ra[#X1"_2"],ra[#X2"_2"]);
#define h3_(X1,X2,X3) hist<1,0,0> \
  h_##X1##_##X2##_##X3(#X1"-"#X2"-"#X3,ra[#X1"_2"],ra[#X2"_2"],ra[#X3"_2"]);

#define hj_(X) auto h_jet_##X = reserve<hist<1>>(njets_expected+1); \
  for (unsigned i=0; i<njmax; ++i) { \
    const auto name = cat("jet",i+1,"_"#X); \
    h_jet_##X.emplace_back(name,ra[name]); \
  }

#define IMPL_HIST_DEFS
#include STR(IMPL)
#undef IMPL_HIST_DEFS

  // ================================================================

  std::vector<fj::PseudoJet> particles;
  Int_t prev_id = -1;

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting \033[36m" << njets_expected
       << "\033[0m or more jets per event\n" << endl;

  Higgs2diphoton Hdecay;
  std::pair<TLorentzVector,TLorentzVector> diphoton;
  boost::optional<TLorentzVector> higgs;

  size_t ncount_total = 0, num_events = 0;

  // LOOP ===========================================================
  using counter = ivanp::timed_counter<Long64_t>;
  for (counter ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    for (unsigned i=_weights.size(); i; ) { // get weights
      --i; bin_t::weights[i] = *_weights[i];
    }

    // Keep track of multi-entry events -----------------------------
    nlo_bin::current_id = *_id;
    const bool new_id = (prev_id != nlo_bin::current_id);
    if (new_id) {
      prev_id = nlo_bin::current_id;
      ncount_total += ( _ncount ? **_ncount : 1);
      ++num_events;
    }
    // --------------------------------------------------------------

    const size_t np = *_nparticle;
    particles.clear();
    higgs = boost::none;

    // Read particles -----------------------------------------------
    unsigned n22 = 0; // number of photons
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

    // const double H_mass = higgs->M();
    const double A1_eta = A1->Eta();
    const double A2_eta = A2->Eta();

    if (A1_pT < 0.35*125.) passes_photon_cuts = false; else
    if (A2_pT < 0.25*125.) passes_photon_cuts = false; else
    if (photon_eta_cut(std::abs(A1_eta))) passes_photon_cuts = false; else
    if (photon_eta_cut(std::abs(A2_eta))) passes_photon_cuts = false;

    cat_bin::id<photon_cuts>() = passes_photon_cuts;

#define IMPL_CAT
#include STR(IMPL)
#undef IMPL_CAT

    h_Njets.fill_bin(njets);

    if (njets < njets_expected) continue; // at least expected number of jets

    // Fill histograms ----------------------------------------------

#define IMPL_HIST_FILL
#include STR(IMPL)
#undef IMPL_HIST_FILL

  } // END EVENT LOOP
  // ================================================================

  cout << "Processed events: " << num_events << endl;
  // TEST(ncount_total)

  // Open output root file for histograms
  auto fout = std::make_unique<TFile>(ofname,"recreate");
  if (fout->IsZombie()) return 1;
  TDirectory *dir = fout.get();

  auto h_Njets_integrated = h_Njets;
  h_Njets_integrated.integrate_left();

#define CATEGORY_TOP(r, data, elem) \
  cat_bin::id<elem>() = 0; \
  for (const char* dir_name : enum_traits<elem>::all_str()) { \
    dir = dir->mkdir(dir_name);

#define CATEGORY_BOT(r, data, elem) \
    dir = dir->GetMotherDir(); \
    ++cat_bin::id<elem>(); \
  }

  // write root historgrams
  wi = 0;
  for (const auto& _w : _weights) {
    dir = dir->mkdir(_w.GetBranchName());
    cout << dir->GetName() << endl;

    BOOST_PP_SEQ_FOR_EACH(CATEGORY_TOP,,CATEGORIES)

      dir->cd();

      using ivanp::root::to_root;
      using ivanp::root::slice_to_root;

      auto* h_Njets_excl = to_root(h_Njets,"Njets_excl");
      excl_labels(h_Njets_excl,true);
      auto* h_Njets_incl = to_root(h_Njets_integrated,"Njets_incl");
      h_Njets_incl->SetEntries( h_Njets_excl->GetEntries() );
      excl_labels(h_Njets_incl,false);

      for (auto& h : hist<1>::all) to_root(*h,h.name);
      for (auto& h : hist<1,0>::all) {
        const auto vars = ivanp::rsplit<1>(h.name,'-');
        slice_to_root(*h,vars[0],vars[1]);
      }

    BOOST_PP_SEQ_FOR_EACH(CATEGORY_BOT,,BOOST_PP_SEQ_REVERSE(CATEGORIES))

    dir = dir->GetMotherDir();
    ++wi;
  }

  fout->cd();
  TH1D *h_N = new TH1D("N","N",1,0,1);
  h_N->SetBinContent(1,ncount_total);
  h_N->SetEntries(num_events);

  fout->Write();

  note("FastJet",jet_def.description_no_recombiner());
  note("Jet cuts",
    cat("pT > ",jet_pt_cut," && ","eta < ",jet_eta_cut));
  note("Photon cuts",
    "(eta < 1.37 || 1.52 < eta) && (eta < 2.37)\n"
    "pT1 > 0.35*mH && pT2 > 0.25*mH");
  note("Input files",
    std::accumulate(std::next(ntuples.begin()),ntuples.end(),
      std::string(ntuples.front()),
      [](std::string s, const auto& x){ return s + '\n' + x; }
    ) );

#define IMPL_INFO
#include STR(IMPL)
#undef IMPL_INFO

  cout << "\n" << tc::green << "Output" << tc::reset
       << ": " << fout->GetName() << endl;

  return 0;
}

