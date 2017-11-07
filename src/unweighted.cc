#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <chrono>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1.h>

#include <fastjet/ClusterSequence.hh>

#include "termcolor.hpp"

#include "float_or_double_reader.hh"
#include "program_options.hh"
#include "timed_counter.hh"
#include "Higgs2diphoton.hh"

#define TEST(var) \
  std::cout << tc::cyan << #var << tc::reset << " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
namespace tc = termcolor;
using namespace ivanp;
namespace fj = fastjet;

std::ostream& operator<<(std::ostream& out, const std::exception& e) {
  return out << tc::red << e.what() << tc::reset;
}

inline bool photon_eta_cut(double abs_eta) noexcept {
  return (1.37 < abs_eta && abs_eta < 1.52) || (2.37 < abs_eta);
}

template <typename... Args>
inline void write(const char* name, Args&&... args) {
  TNamed(name,cat(args...).c_str()).Write();
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char *ofname = nullptr, *tree_name = "t3";
  double max_weight_sampling = 1e3;
  double jet_pt_cut = 30., jet_eta_cut = 4.4;
  unsigned njets_expected;
  bool no_photon_cuts = false;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT ntuples",req(),pos())
      (ofname,'o',"output file name")
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (max_weight_sampling,{"-w","--max-weight"},
       cat("maximum sampling weight [",max_weight_sampling,']'))
      (njets_expected,'j',"expected number of jets",req())
      (jet_pt_cut,"--jet-pt-cut",cat('[',jet_pt_cut,']'))
      (jet_eta_cut,"--jet-eta-cut",cat('[',jet_eta_cut,']'))
      (no_photon_cuts,"--no-photon-cuts")
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  auto rand = [=]{ // mersenne twister random number generator
    static std::mt19937 gen(
      std::chrono::system_clock::now().time_since_epoch().count() );
    static std::uniform_real_distribution<double> dist(0.,max_weight_sampling);
    return dist(gen);
  };

  TChain chain(tree_name);
  cout << tc::cyan << "Input ntuples" << tc::reset << ':' << endl;
  for (const char* name : ifnames) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);
  TTreeReaderValue<Double_t> _weight(reader,"weight2");

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  // Open output file
  TFile *fout;
  TTree *tree;
  TH1D *h_weight;

  unsigned char np;
  // int pid[2];
  double p[4][5];

  if (ofname) {
    fout = new TFile(ofname,"recreate","unweighted events",109);
    if (fout->IsZombie()) return 1;
    cout << tc::cyan << "Output file" << tc::reset << ": "
         << fout->GetName() << endl;

    h_weight = new TH1D("weight","weight",1000,0,std::log10(max_weight_sampling));

    tree = new TTree("events","events");
    tree->Branch("np",&np);
    // tree->Branch("pid",pid,"pid[np]/I");
    tree->Branch("px",p[0], "px[np]/D");
    tree->Branch("py",p[1], "py[np]/D");
    tree->Branch("pz",p[2], "pz[np]/D");
    tree->Branch( "E",p[3],  "E[np]/D");
  }

  const fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() <<'\n'<< endl;
  std::vector<fj::PseudoJet> particles;

  Higgs2diphoton Hdecay;

  // LOOP ===========================================================
  unsigned long events_selected = 0, events_used = 0;
  double max_weight_in_file = 0.,
         max_weight_jet_cuts = 0,
         max_weight_photon_cuts = 0;
  using cnt = timed_counter<Long64_t>;
  for (cnt ent(events_used=reader.GetEntries(true)); reader.Next(); ++ent) {
    np = *_nparticle;
    particles.clear();
    for (decltype(np) i=0; i<np; ++i) {
      if (_kf[i]==25) {
        p[0][0] = _px[i];
        p[1][0] = _py[i];
        p[2][0] = _pz[i];
        p[3][0] = _E [i];
      } else {
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
      }
    }

    const double weight = *_weight;
    if (weight > max_weight_in_file) max_weight_in_file = weight;

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(jet_pt_cut); // apply pT cut

    // apply eta cut
    for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
      if (std::abs(it->eta()) > jet_eta_cut) fj_jets.erase(it);
      else ++it;
    }

    // resulting number of jets
    const unsigned njets = fj_jets.size();
    if (njets < njets_expected) continue;

    // sort by pT
    std::sort( fj_jets.begin(), fj_jets.end(),
      [](const fj::PseudoJet& a, const fj::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    // --------------------------------------------------------------

    if (weight > max_weight_jet_cuts) max_weight_jet_cuts = weight;

    if (!ofname) continue;

    // Photons ------------------------------------------------------
    if (!no_photon_cuts) {
      const TLorentzVector Higgs(p[0][0],p[1][0],p[2][0],p[3][0]);
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

      if (weight > max_weight_photon_cuts) max_weight_photon_cuts = weight;
    }
    // --------------------------------------------------------------

    h_weight->Fill(std::log10(weight));

    if (rand() > weight) continue;
    ++events_selected;

    for (unsigned j=0; j<njets; ++j) {
      p[0][j+1] = fj_jets[j][0];
      p[1][j+1] = fj_jets[j][1];
      p[2][j+1] = fj_jets[j][2];
      p[3][j+1] = fj_jets[j][3];
    }

    tree->Fill();
  } // end loop
  cout << endl;
  // ================================================================

  TEST( max_weight_sampling )
  TEST( max_weight_in_file )
  TEST( max_weight_jet_cuts )
  TEST( max_weight_photon_cuts )
  TEST( events_used )

  if (ofname) {
    TEST( events_selected )

    const double efficiency = 100.*double(events_selected)/double(events_used);
    cout << tc::cyan << "efficiency" << tc::reset << " = "
         << efficiency << '%' << std::endl;

    write("max_weight_sampling",max_weight_sampling);
    write("max_weight_in_file",max_weight_in_file);
    write("max_weight_jet_cuts",max_weight_jet_cuts);
    write("max_weight_photon_cuts",max_weight_photon_cuts);
    write("events_used",events_used);
    write("events_selected",events_selected);
    write("efficiency",std::setprecision(3),efficiency,'%');

    fout->Write(0,TObject::kOverwrite);
    delete fout;
  }
}
