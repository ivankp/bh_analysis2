// Written by Ivan Pogrebnyak

#include <iostream>
#include <cstring>
#include <algorithm>
#include <memory>

#include <boost/optional.hpp>

#include <TFile.h>
// #include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "math.hh"
#include "timed_counter.hh"
#include "catstr.hh"
#include "exception.hh"
#include "float_or_double_reader.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;

using ivanp::cat;

using namespace ivanp::math;

int main(int argc, char* argv[]) {
  // Open input ntuples root file ===================================
  TChain chain("BHSntuples");
  cout << "\033[36mInput ntuples\033[0m:" << endl;
  for (int i=1; i<argc; ++i) {
    if (!chain.Add(argv[i],0)) return 1;
    cout << "  " << argv[i] << endl;
  }
  cout << endl;

  // Set up branches for reading
  TTreeReader reader(&chain);

  TTreeReaderValue<Int_t> _id(reader,"id");
  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  boost::optional<TTreeReaderValue<Char_t>> _part;
  for ( auto bo : *reader.GetTree()->GetListOfBranches() ) {
    if (!strcmp(bo->GetName(),"part")) _part.emplace(reader,"part");
  }

  float_or_double_value_reader _weight(reader,"weight");

  // LOOP ===========================================================
  using tc = ivanp::timed_counter<Long64_t>;
  for (tc ent(reader.GetEntries(true)); reader.Next() && ent < 10; ++ent) {
    test( *_id )
    test( *_weight )
  }

  return 0;
}

