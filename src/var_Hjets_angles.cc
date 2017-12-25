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
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "catstr.hh"
// #include "exception.hh"
#include "float_or_double_reader.hh"
#include "Higgs2diphoton.hh"
// #include "comprehension.hh"
#include "program_options.hh"
#include "parse_args/jetdef.hh"
#include "tc_msg.hh"
// #include "string_alg.hh"

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

constexpr bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname = nullptr;
  const char* tree_name = "t3";
  const char* weight_branch_name = "weight";
  fj::JetDefinition jet_def;
  unsigned njets_expected;
  double jet_pt_cut = 30., jet_eta_cut = 4.4;
  bool no_photon_cuts = false;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT BH ntuples",req(),pos())
      (ofname,'o',"output file name",req())
      (tree_name,{"-t","--tree"},cat("input TTree name [",tree_name,']'))
      (weight_branch_name,{"-w","--weight"},
       cat("weight branch [",weight_branch_name,']'))
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

  // Open input ntuples root file ===================================
  TChain chain(tree_name);
  info("Input ntuples");
  for (const char* name : ifnames) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);

  // TTreeReaderValue<Int_t> _id(reader,"id");
  boost::optional<TTreeReaderValue<Int_t>> _nparticle;
  boost::optional<TTreeReaderArray<Int_t>> _kf;

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  float_or_double_value_reader _weight(reader,weight_branch_name);

  boost::optional<TTreeReaderValue<UChar_t>> _ncount;

#define OPT_BRANCH(NAME) \
  if (!strcmp(bo->GetName(),#NAME)) _##NAME.emplace(reader,#NAME);

  for (auto bo : *reader.GetTree()->GetListOfBranches()) {
    OPT_BRANCH(ncount)
    OPT_BRANCH(nparticle)
    OPT_BRANCH(kf)
    if (_ncount && _nparticle && _kf) break;
  }

#undef OPT_BRANCH

  // Output file ====================================================
  TFile fout(ofname,"recreate","",109);
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TTree *tout = new TTree("angles","");
  Int_t ncount;
  double weight, hj_mass, cos_theta;
  tout->Branch("weight",&weight);
  if (_ncount) tout->Branch("ncount",&ncount);
  tout->Branch("hj_mass",&hj_mass);
  tout->Branch("cos_theta",&cos_theta);

  if (_kf) {
    write("FastJet",jet_def.description());
    write("Jet cuts","pT > ",jet_pt_cut," GeV\neta < ",jet_eta_cut);
    write("Photon cuts","pT_1 > 0.35 M_Higgs\npT_2 > 0.25 M_Higgs\n"
                        "eta < 1.37 or 1.52 < eta < 2.37");
  }
  // TODO: else read from input file

  // ================================================================
  std::vector<fj::PseudoJet> particles;
  particles.reserve(njets_expected+1);

  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  cout << "Expecting " << njets_expected
       << " or more jets per event\n" << endl;

  TLorentzVector Higgs, jet1;
  Higgs2diphoton Hdecay;

  // LOOP ===========================================================
  using cnt = ivanp::timed_counter<Long64_t>;
  for (cnt ent(reader.GetEntries(true)); reader.Next(); ++ent) {

    if (_kf) {
      // Read particles -----------------------------------------------
      const unsigned np = **_nparticle;
      particles.clear();
      for (unsigned i=0; i<np; ++i) {
        if ((*_kf)[i]==25) {
          Higgs.SetPxPyPzE(_px[i],_py[i],_pz[i],_E[i]);
        } else {
          particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
        }
      }
      // --------------------------------------------------------------

      // Cluster jets -----------------------------------------------
      auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
        .inclusive_jets(jet_pt_cut); // apply pT cut
      // apply eta cut
      for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
        if (std::abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
        else ++it;
      }
      // resulting number of jets
      const unsigned njets = fj_jets.size();
      // require expected number of jets
      if (njets < njets_expected) continue;
      // sort by pT
      std::sort( fj_jets.begin(), fj_jets.end(),
        [](const fj::PseudoJet& a, const fj::PseudoJet& b){
          return ( a.pt() > b.pt() );
        });

      // Photons ----------------------------------------------------
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
      // ------------------------------------------------------------

      jet1.SetPxPyPzE(
        fj_jets.front().px(),
        fj_jets.front().py(),
        fj_jets.front().pz(),
        fj_jets.front().e()
      );

    } else {
      Higgs.SetPxPyPzE(_px[0],_py[0],_pz[0],_E[0]);
      jet1 .SetPxPyPzE(_px[1],_py[1],_pz[1],_E[1]);
    }

    // --------------------------------------------------------------
    weight = *_weight;
    if (_ncount) ncount = **_ncount;
    // --------------------------------------------------------------

    const auto Q = Higgs + jet1;
    const double Q2 = Q*Q;
    hj_mass = sqrt(Q2);

    const TLorentzVector Z(0,0,Q.E(),Q.Pz());
    const auto ell = ((Q*jet1)/Q2)*Higgs - ((Q*Higgs)/Q2)*jet1;

    cos_theta = (ell*Z) / sqrt(sq(ell)*sq(Z));

    tout->Fill();
  } // END EVENT LOOP
  // ================================================================
  cout.imbue(std::locale(""));
  info("Events selected", tout->GetEntries());

  fout.Write(0,TObject::kOverwrite);

  return 0;
}
