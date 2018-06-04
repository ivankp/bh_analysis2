// Written by Ivan Pogrebnyak

#include <iostream>
#include <fstream>
#include <algorithm>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <TFile.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "float_or_double_reader.hh"
#include "timed_counter.hh"
#include "catstr.hh"
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

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname = nullptr;
  const char* tree_name = "t3";
  unsigned prec = 15;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT ntuples",req(),pos())
      (ofname,'o',"output file name",req())
      (tree_name,{"-t","--tree"},cat("input TTree name [",tree_name,']'))
      (prec,'p',cat("precision [",prec,']'))
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

  TTreeReader reader(&chain);
  TTreeReaderValue<Int_t> nparticle(reader,"nparticle");
  float_or_double_value_reader weight(reader,"weight2");
  float_or_double_array_reader px(reader,"px");
  float_or_double_array_reader py(reader,"py");
  float_or_double_array_reader pz(reader,"pz");
  float_or_double_array_reader E (reader,"E" );

  // Output file ====================================================
  std::ofstream outf(ofname, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  out.push(boost::iostreams::bzip2_compressor());
  out.push(outf);
  out << std::scientific << std::setprecision(prec);

  // LOOP ===========================================================
  using counter = ivanp::timed_counter<Long64_t>;
  for (counter ent(reader.GetEntries(true)); reader.Next(); ++ent) {

    out << *weight;
    for (auto n=*nparticle, i=0; i<n; ++i)
      out << ',' << E[i] << ',' << px[i] << ',' << py[i] << ',' << pz[i];
    out << '\n';

  } // END EVENT LOOP
  // ================================================================

  return 0;
}
