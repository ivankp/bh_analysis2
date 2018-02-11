#include <iostream>
#include <string>
#include <map>
#include <set>

#include <TFile.h>
#include <TDirectory.h>
#include <TH1.h>

#include "tkey.hh"
#include "catstr.hh"

using std::cout;
using std::endl;
using std::cerr;

template <typename T>
inline T* get(TDirectory* dir, const char* name) {
  T *obj = nullptr;
  dir->GetObject(name,obj);
  if (obj) return obj;
  else throw std::runtime_error(
    ivanp::cat("No object ",name," in ",dir->GetName()));
}

const char* current_file_name;
std::map<TH1*,std::set<const char*>> hist_tracker;
double N_scale;

void loop(TDirectory* dir) { // LOOP
  for (TKey& key : get_keys(dir)) {
    const TClass* key_class = get_class(key);
    if (inherits_from<TH1>(key_class)) { // HIST
      TH1 *in  = key_cast<TH1>(key),
          *out = nullptr;
      gDirectory->GetObject(in->GetName(),out);
      if (out) {
        const auto ent = out->GetEntries() + in->GetEntries();
        out->Add(in,N_scale);
        out->SetEntries(ent);
      } else {
        in->Scale(N_scale);
        out = static_cast<TH1*>(in->Clone());
      }
      // keep track of histograms
      hist_tracker[out].insert(current_file_name);

    } else if (inherits_from<TDirectory>(key_class)) { // DIR
      TDirectory *dir = key_cast<TDirectory>(key),
                 *out = gDirectory;
      const char* name = dir->GetName();
      TDirectory *out2 = out->GetDirectory(name);
      if (!out2) out2 = out->mkdir(name);
      out2->cd();
      loop(dir);
      out->cd();
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc<3) {
    cout << "usage: " << argv[0] << " NLO.root B.root RS.root ..." << endl;
    return 1;
  }

  // Output file
  TFile fout(argv[1],"recreate");
  if (fout.IsZombie()) return 1;
  cout << "Output file: " << fout.GetName() << endl;

  for (int fi=2; fi<argc; ++fi) {
    TFile f((current_file_name = argv[fi]));
    if (f.IsZombie()) return 1;
    cout << "Input file: " << f.GetName() << endl;

    // Number of events
    const double N = get<TH1D>(&f,"N")->GetAt(1);
    N_scale = 1./N;
    cout << "Events: " << N << endl;

    fout.cd();
    loop(&f); // LOOP
  }

  // Set N to 1
  get<TH1D>(&fout,"N")->SetAt(1,1);

  // Write and close output root file
  fout.Write();
  cout << "\n\033[32mWrote file\033[0m: " << fout.GetName() << endl;

  std::map<std::set<const char*>,std::vector<TH1*>> hist_tracker2;
  for (const auto& h : hist_tracker)
    hist_tracker2[h.second].push_back(h.first);
  if (hist_tracker2.size()>1) {
    cerr << "\n\033[33m"
      "Input files contained different sets of histograms"
      "\033[0m" << endl;
  }
}
