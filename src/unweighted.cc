#include <iostream>
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

#include "termcolor.hpp"

#include "float_or_double_reader.hh"
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
  const char *ofname, *tree_name = "t3";
  double max_weight = 1e3;

  try {
    using namespace ivanp::po;
    if (program_options()
      (ifnames,'i',"input ROOT ntuples",req(),pos())
      (ofname,'o',"output file name",req())
      (tree_name,"--tree",cat("input TTree name [",tree_name,']'))
      (max_weight,{"-w","--max-weight"},cat("maximum weight [",max_weight,']'))
      .parse(argc,argv,true)) return 0;
  } catch (const std::exception& e) {
    cerr << e << endl;
    return 1;
  }

  auto rand = [=]{ // mersenne twister random number generator
    static std::mt19937 gen(
      std::chrono::system_clock::now().time_since_epoch().count() );
    static std::uniform_real_distribution<double> dist(0.,max_weight);
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

  TTreeReaderValue<Int_t> _nparticle(reader,"nparticle");
  TTreeReaderArray<Int_t> _kf(reader,"kf");
  TTreeReaderValue<Double_t> _weight(reader,"weight2");

  float_or_double_array_reader _px(reader,"px");
  float_or_double_array_reader _py(reader,"py");
  float_or_double_array_reader _pz(reader,"pz");
  float_or_double_array_reader _E (reader,"E" );

  // Open output file
  TFile fout(ofname,"recreate");
  if (fout.IsZombie()) return 1;
  cout << tc::cyan << "Output file" << tc::reset << ':' << endl;

  auto* h_weight = new TH1D("weight","weight",1000,0,std::log10(max_weight));
  max_weight = 0.;

  unsigned char np;
  int pid[5];
  double p[4][5];

  TTree *tree = new TTree("events","events");
  tree->Branch("np",&np);
  tree->Branch("pid",pid,"pid[np]/I");
  tree->Branch("px",p[0],"px[np]/D");
  tree->Branch("py",p[1],"py[np]/D");
  tree->Branch("pz",p[2],"pz[np]/D");
  tree->Branch( "E",p[3], "E[np]/D");

  // LOOP ===========================================================
  unsigned long selected = 0, total;
  for (timed_counter<Long64_t> ent(total=reader.GetEntries(true)); reader.Next(); ++ent) {
    const auto weight = *_weight;
    h_weight->Fill(std::log10(weight));
    if (weight > max_weight) max_weight = weight;

    if (rand() > weight) continue;
    ++selected;

    np = *_nparticle;
    for (decltype(np) i=0; i<np; ++i) {
      pid [i] = _kf[i];
      p[0][i] = _px[i];
      p[1][i] = _py[i];
      p[2][i] = _pz[i];
      p[3][i] = _E [i];
    }

    tree->Fill();
  }
  cout << endl;

  TEST( max_weight )
  TEST( total )
  TEST( selected )

  cout << tc::cyan << "efficiency" << tc::reset << ": "
       << ( 100.*double(selected)/double(total) ) << std::endl;

  fout.Write(0,TObject::kOverwrite);
}
