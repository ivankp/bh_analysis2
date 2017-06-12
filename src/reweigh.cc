// Written by Ivan Pogrebnyak

#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include <boost/program_options.hpp>

#include <TFile.h>
#include <TTree.h>
#include <TNamed.h>

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
               "recreate","weights",109);
    {
      TNamed("central scale",scales_pretty.at(scale_name).c_str()).Write();
      TNamed("PDF",pdf_name.c_str()).Write();
      TNamed("ntuple",ntuple.c_str()).Write();
      const time_t current_time = std::time(nullptr);
      TNamed("time stamp",std::ctime(&current_time)).Write();
    }
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
    for (const auto& i : combs) {
      sd.scales.emplace_back(i.first,i.second);
    }

    reweighter rew(*tin,sd,pdf_name,do_pdf_variations);

    for (unsigned i=0; i<combs.size(); ++i) {
      tout->Branch(
        cat( "MUF",fracs[sd.scales_fac[combs[i].first]],
            "_MUR",fracs[sd.scales_ren[combs[i].second]],
            "_PDF",rew.pdf_id(0) ).c_str(),
        &rew[i]);
    }
    if (do_pdf_variations) {
      for (unsigned i=combs.size(); i<rew.size(); ++i) {
        tout->Branch(
          cat( "MUF1_MUR1_PDF", rew.pdf_id(i+1-combs.size()) ).c_str(),
          &rew[i]);
      }
    }

    // LOOP ===========================================================
    using tc = ivanp::timed_counter<Long64_t>;
    for (tc ent(tin->GetEntries()); !!ent; ++ent) {
      tin->GetEntry(ent);

      rew();

      tout->Fill();
    }

    tout->Write("",TObject::kOverwrite);
    cout << "\n\033[32mWrote\033[0m: " << fout.GetName() << endl;
    cout << "Compression factor: " << fout.GetCompressionFactor() << endl;
  }

  return 0;
}

