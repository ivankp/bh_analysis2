#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include "termcolor.hpp"

#include "program_options.hh"
#include "timed_counter.hh"

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
  const char *ofname, *tree_name = "events";
  unsigned nevents = 10000;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT files",req(),pos())
      (ofname,'o',"output file name",req())
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (nevents,'n',cat("number of events [",nevents,']'))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  TChain chain(tree_name);
  cout << tc::cyan << "Input files" << tc::reset << ':' << endl;
  for (const char* name : ifnames) {
    if (!chain.Add(name,0)) return 1;
    cout << "  " << name << endl;
  }
  cout << endl;

  // Set up branches for reading
  UChar_t np;
  Int_t pid[5];
  Double_t p[4][5];
  
  chain.SetBranchAddress("np",&np);
  chain.SetBranchAddress("pid",pid);
  chain.SetBranchAddress("px",p[0]);
  chain.SetBranchAddress("py",p[1]);
  chain.SetBranchAddress("pz",p[2]);
  chain.SetBranchAddress( "E",p[3]);

  // Open output file
  std::ofstream f(ofname);

  // LOOP ===========================================================
  for (timed_counter<Long64_t> ent(nevents); !!ent; ++ent) {
    chain.GetEntry(ent);
    for (UChar_t i=0; i<np; ++i) {
      if (i) f << ' ';
      f << pid[i];
      for (unsigned j=0; j<4; ++j)
        f << ' ' << p[j][i];
    }
    f << endl;
  }
}

