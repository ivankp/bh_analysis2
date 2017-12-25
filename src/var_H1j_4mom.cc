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

#include "math.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "float_or_double_reader.hh"
#include "Higgs2diphoton.hh"
#include "program_options.hh"
#include "tc_msg.hh"

#define TEST(VAR) \
  std::cout << tc::cyan << #VAR << tc::reset << " = " << VAR << std::endl;

#define _STR(S) #S
#define STR(S) _STR(S)

using std::cout;
using std::cerr;
using std::endl;
namespace tc = termcolor;
using namespace ivanp;
using namespace ivanp::math;

struct error : std::runtime_error {
  using std::runtime_error::runtime_error;
  template <typename... Args>
  error(const Args&... args): error(cat(args...)) { }
};

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
  double jet_pt_cut = 30., jet_eta_cut = 4.4;
  bool no_photon_cuts = false;
  const char* weight_branch_name = "weight";

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT BH ntuples",req(),pos())
      (ofname,'o',"output file name",req())
      (tree_name,{"-t","--tree"},cat("input TTree name [",tree_name,']'))
      (jet_pt_cut,"--jet-pt-cut",cat('[',jet_pt_cut,"] GeV"))
      (jet_eta_cut,"--jet-eta-cut",cat('[',jet_eta_cut,']'))
      (no_photon_cuts,"--no-photon-cuts")
      (weight_branch_name,{"-w","--weight"},
       cat("weight branch [",weight_branch_name,']'))
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
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  float_or_double_value_reader _weight(reader,weight_branch_name);

  boost::optional<TTreeReaderValue<Int_t>> _ncount;
  for (auto bo : *reader.GetTree()->GetListOfBranches()) {
    if (!strcmp(bo->GetName(),"ncount")) _ncount.emplace(reader,"ncount");
    if (_ncount) break;
  }

  // Output file ====================================================
  TFile fout(ofname,"recreate","",109);
  info("Output file",fout.GetName());
  if (fout.IsZombie()) return 1;
  fout.cd();

  TTree *tout = new TTree("Hj","");
  UChar_t ncount;
  Double_t weight,
           px[2], py[2], pz[2], E[2];
  tout->Branch("weight",&weight);
  if (_ncount) tout->Branch("ncount",&ncount);
  tout->Branch("px",px,"px[2]/D");
  tout->Branch("py",py,"py[2]/D");
  tout->Branch("pz",pz,"pz[2]/D");
  tout->Branch("E" ,E , "E[2]/D");

  write("Jet cuts","pT > ",jet_pt_cut," GeV\neta < ",jet_eta_cut);
  if (!no_photon_cuts) {
    write("Photon cuts","pT_1 > 0.35 M_Higgs\npT_2 > 0.25 M_Higgs\n"
                        "eta < 1.37 or 1.52 < eta < 2.37");
  }

  // ================================================================
  TLorentzVector Higgs, jet1;
  Higgs2diphoton Hdecay;

  // LOOP ===========================================================
  using cnt = ivanp::timed_counter<Long64_t>;
  for (cnt ent(reader.GetEntries(true)); reader.Next(); ++ent) {

    // Read particles -----------------------------------------------
    const int np = *_nparticle;
    if (np!=2)
      throw error(np," particles in event ",ent);
    if (_kf[0]!=25)
      throw error("Higgs is not first particle in event ",ent);

    Higgs.SetPxPyPzE(_px[0],_py[0],_pz[0],_E[0]);
    jet1 .SetPxPyPzE(_px[1],_py[1],_pz[1],_E[1]);
    // --------------------------------------------------------------

    if (jet1.Pt() < jet_pt_cut) continue;
    if (jet1.Eta() > jet_eta_cut) continue;

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

    // --------------------------------------------------------------
    weight = *_weight;
    if (_ncount) ncount = **_ncount;

    for (int i=0; i<np; ++i) {
      px[i] = _px[i];
      py[i] = _py[i];
      pz[i] = _pz[i];
      E [i] = _E [i];
    }
    // --------------------------------------------------------------

    tout->Fill();
  } // END EVENT LOOP
  // ================================================================
  cout.imbue(std::locale(""));
  info("Events selected", tout->GetEntries());

  fout.Write(0,TObject::kOverwrite);

  return 0;
}
