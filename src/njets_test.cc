#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <memory>
#include <array>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TLorentzVector.h>

#include <fastjet/ClusterSequence.hh>

#include "math.hh"
#include "timed_counter.hh"
#include "float_or_double_reader.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using ivanp::math::sq;

namespace fj = fastjet;

int main(int argc, char* argv[]) {
  auto f = std::make_unique<TFile>(argv[1],"read");
  if (f->IsZombie()) return 1;

  TTreeReader reader("t3",f.get());

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );
  
  float_or_double_value_reader _weight(reader, "weight");

  std::array<double,7> npar{}, npar_30{}, njets{},
                       epar{}, epar_30{}, ejets{};

  const fj::JetDefinition jet_def(fj::antikt_algorithm,0.4);
  fastjet::ClusterSequence::print_banner(); // get it out of the way
  cout << jet_def.description() << endl;
  std::vector<fj::PseudoJet> particles;

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    const double weight = *_weight; // Read weight

    // Read particles -----------------------------------------------
    particles.clear();
    particles.reserve(*_nparticle);
    for (size_t i=0, n=_kf.GetSize(); i<n; ++i)
      if (_kf[i] != 25)
        particles.emplace_back(_px[i],_py[i],_pz[i],_E[i]);
    // --------------------------------------------------------------
    const size_t np = particles.size();

    npar[np] += weight;
    ++epar[np];
    
    int np_30 = 0;
    for (const auto& p : particles)
      if ( p.pt() > 30. ) ++np_30;
    npar_30[np_30] += weight;
    ++epar_30[np_30];

    // Cluster jets -------------------------------------------------
    auto fj_jets = fj::ClusterSequence(particles,jet_def) // cluster
      .inclusive_jets(30.); // apply pT cut
    // apply eta cut
    for (auto it=fj_jets.begin(); it!=fj_jets.end(); ) {
      if (abs(it->eta()) > 4.4) fj_jets.erase(it);
      else ++it;
    }
    // sort by pT
    std::sort( fj_jets.begin(), fj_jets.end(),
      [](const fj::PseudoJet& a, const fj::PseudoJet& b){
        return ( a.pt() > b.pt() );
      });
    // resulting number of jets
    const unsigned nj = fj_jets.size();
    // --------------------------------------------------------------

    njets[nj] += weight;
    ++ejets[nj];
  }

  cout.precision(3);
  for (unsigned i=0; i<njets.size(); ++i) {
    cout << i << ": "
         << setw(10) << epar[i] << ' '
         << setw(10) << epar_30[i] << " "
         << setw(10) << ejets[i] << endl;
  }
  cout << endl;
  for (unsigned i=0; i<njets.size(); ++i) {
    cout << i << ": "
         << setw(10) << npar[i] << ' '
         << setw(10) << npar_30[i] << " "
         << setw(10) << njets[i] << endl;
  }

  return 0;
}
