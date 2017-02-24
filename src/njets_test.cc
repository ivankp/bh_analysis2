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

  std::array<double,7> npar{}, njets{}, njets_30{},
                       epar{}, ejets{}, ejets_30{};
  std::vector<std::array<double,4>> jets;

  // LOOP ***********************************************************
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next(); ++ent) {
    const double weight = *_weight; // Read weight

    const int np = *_nparticle;
    // Read particles -----------------------------------------------
    jets.clear();
    for (int i=0; i<np; ++i)
      if (_kf[i] != 25)
        jets.push_back({_px[i],_py[i],_pz[i],_E[i]});
    // --------------------------------------------------------------

    npar[np] += weight;
    ++epar[np];
    njets[jets.size()] += weight;
    ++ejets[jets.size()];
    
    int nj_30 = 0;
    for (const auto& j : jets)
      if ( std::sqrt(sq(j[0],j[1])) > 30. ) ++nj_30;
    njets_30[nj_30] += weight;
    ++ejets_30[nj_30];

  }

  cout.precision(3);
  for (unsigned i=0; i<njets.size(); ++i) {
    cout << i << ": "
         << setw(10) << epar[i] << ' '
         << setw(10) << ejets[i] << " "
         << setw(10) << ejets_30[i] << endl;
  }
  cout << endl;
  for (unsigned i=0; i<njets.size(); ++i) {
    cout << i << ": "
         << setw(10) << npar[i] << ' '
         << setw(10) << njets[i] << " "
         << setw(10) << njets_30[i] << endl;
  }

  return 0;
}
