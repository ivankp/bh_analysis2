#include <iostream>
#include <iomanip>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "termcolor.hpp"

#include "program_options.hh"
// #include "timed_counter.hh"

#define TEST(var) \
  std::cout << tc::cyan << #var << tc::reset << " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
namespace tc = termcolor;
using namespace ivanp;

std::ostream& operator<<(std::ostream& out, const std::exception& e) {
  return out << tc::red << e.what() << tc::reset;
}

int main(int argc, char* argv[]) {
  std::vector<const char*> ifnames;
  const char *ofname = nullptr, *tree_name = "events";
  unsigned nevents = 0;
  unsigned prec = 15;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT files",req(),pos())
      (ofname,'o',"output file name")
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (nevents,'n',cat("number of events [",nevents,']'))
      (prec,'p',cat("precision [",prec,']'))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TChain chain(tree_name);
  // cout << tc::cyan << "Input files" << tc::reset << ':' << endl;
  for (const char* name : ifnames) {
    if (!chain.Add(name,0)) return 1;
    // cout << "  " << name << endl;
  }
  // cout << endl;

  // Set up branches for reading
  UChar_t np;
  // Int_t pid[5];
  Double_t p[4][5];
  
  chain.SetBranchAddress("np",&np);
  // chain.SetBranchAddress("pid",pid);
  chain.SetBranchAddress("px",p[0]);
  chain.SetBranchAddress("py",p[1]);
  chain.SetBranchAddress("pz",p[2]);
  chain.SetBranchAddress( "E",p[3]);

  // Output stream
  std::ostream& out = ofname ? *new std::ofstream(ofname) : cout;
  out.precision(prec);

  // LOOP ===========================================================
  // for (timed_counter<Long64_t> ent(nevents); !!ent; ++ent) {
  if (!nevents) nevents = chain.GetEntries();
  for (long ent=0; ent<nevents; ++ent) {
    chain.GetEntry(ent);
    for (UChar_t i=0; i<np; ++i) {
      if (i) out << ' ';
      out << p[0][i] << ' '
          << p[1][i] << ' '
          << p[2][i] << ' '
          << p[3][i];
    }
    out << endl;
  }

  if (ofname) delete &out;
}

