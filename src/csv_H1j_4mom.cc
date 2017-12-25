// Written by Ivan Pogrebnyak

#include <iostream>
#include <fstream>
#include <algorithm>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

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

struct error : std::runtime_error {
  using std::runtime_error::runtime_error;
  template <typename... Args>
  error(const Args&... args): error(cat(args...)) { }
};

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char* ofname = nullptr;
  const char* tree_name = "Hj";
  unsigned prec = 15;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT BH ntuples",req(),pos())
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

  // UChar_t ncount;
  Double_t weight,
           px[2], py[2], pz[2], E[2];

  chain.SetBranchAddress("weight",&weight);
  // chain.SetBranchAddress("ncount",&ncount);
  chain.SetBranchAddress("px",px);
  chain.SetBranchAddress("py",py);
  chain.SetBranchAddress("pz",pz);
  chain.SetBranchAddress("E" ,E );

  // Output file ====================================================
  std::ofstream outf(ofname, std::ios_base::out | std::ios_base::binary);
  boost::iostreams::filtering_ostream out;
  out.push(boost::iostreams::bzip2_compressor());
  out.push(outf);
  out << std::scientific << std::setprecision(prec);

  // LOOP ===========================================================
  using cnt = ivanp::timed_counter<Long64_t>;
  for (cnt ent(chain.GetEntries()); !!ent; ++ent) {
    chain.GetEntry(ent);

    out << weight << ','
        // << unsigned(ncount) << ','
        << E[0] << ',' << px[0] << ',' << py[0] << ',' << pz[0] << ','
        << E[1] << ',' << px[1] << ',' << py[1] << ',' << pz[1] << '\n';

  } // END EVENT LOOP
  // ================================================================

  return 0;
}
