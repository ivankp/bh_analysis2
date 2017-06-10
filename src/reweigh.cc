// Written by Ivan Pogrebnyak

#include <iostream>
#include <string>
#include <vector>
// #include <stdexcept>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>

#include "math.hh"
#include "timed_counter.hh"
#include "reweighter.hh"
#include "catstr.hh"

#define test(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
namespace po = boost::program_options;
using ivanp::cat;
using namespace ivanp::math;
#include "scales.hh"

int main(int argc, char* argv[]) {
  // START OPTIONS **************************************************
  std::vector<std::string> ntuples;
  std::string scale_name, pdf_name, output_dir;
  bool do_pdf_variations;

  try {
    // General Options ------------------------------------
    po::options_description desc("Options");
    desc.add_options()
      ("help,h", "produce help message")
      ("input,i", po::value(&ntuples)->multitoken()->required(),
       "*add input ntuples")
      ("output,o", po::value(&output_dir)->required(),
       "*output directory")
      ("scale,s", po::value(&scale_name)->required(),
       "*central scale name")
      ("pdf,p", po::value(&pdf_name)->required(),
       "*PDF set name")
      ("pdf-variations,v", po::bool_switch(&do_pdf_variations),
       "compute weights for PDF uncertainties")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (argc == 1 || vm.count("help")) {
      cout << desc << endl;
      return 0;
    }
    po::notify(vm);
  } catch(std::exception& e) {
    cerr << "\033[31mError: " <<  e.what() <<"\033[0m"<< endl;
    return 1;
  }
  // END OPTIONS ****************************************************

  for (const std::string& ntuple : ntuples) {
    cout << "\033[36mReweighting\033[0m " << ntuple << endl;

    TFile fin(ntuple.c_str());
    TTree *tin = static_cast<TTree*>(fin.Get("t3"));

    TFile fout(cat(output_dir,'/',ntuple.substr(ntuple.rfind('/')+1)).c_str(),
               "recreate");
    TTree *tout = new TTree("weights","");

    scale_defs sd;

    std::vector<double> fracs {1.,2.,0.5};
    for (double x : fracs)
      sd.scale_fcns.emplace_back(
        [x,scale=scales.at(scale_name)](const entry& e){ return x*scale(e); });

    sd.scales_fac = {0,1,2};
    sd.scales_ren = {0,1,2};

    const std::vector<std::pair<unsigned,unsigned>> combs {
      {0,0},{0,1},{1,0},{1,1},{0,2},{2,0},{2,2} };
    sd.scales.reserve(combs.size());
    for (const std::pair<unsigned,unsigned>& i : combs) {
      sd.scales.emplace_back(i.first,i.second);
      tout->Branch(
        cat( "MUF",fracs[i.first],"_MUR",fracs[i.second],"_PDF"/*NUMBER*/).c_str(),
        &sd.scales.back());
    }

    reweighter rew(*tin,pdf_name,sd);

    // LOOP ===========================================================
    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(tin->GetEntries()); !!ent; ++ent) {
      tin->GetEntry(ent);

      rew();

      // test( rew[0] )
      for (unsigned i=0, n=sd.scales.size(); i<n; ++i)
        cout << rew[i] << endl;
    }
  }

  return 0;
}

