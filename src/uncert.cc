#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <stdexcept>

#include <TNamed.h>
#include <TKey.h>
#include <TClass.h>
#include <TFile.h>
#include <TH1.h>

#include <LHAPDF/LHAPDF.h>

#include "utility.hh"
#include "catstr.hh"

#define TEST(var) \
  std::cout <<"\033[36m"<< #var <<"\033[0m"<< " = " << var << std::endl;

using std::cout;
using std::cerr;
using std::endl;
using ivanp::cat;

template <typename...> struct bad_type;

using dirs_t = std::vector<ivanp::named_ptr<TDirectory>>;
using hists_t = std::vector<std::pair< TH1*,
                                       std::vector<std::vector<double>> >>;

auto get_values(const dirs_t& dirs) {
  hists_t hists;
  for (const auto& dir : dirs) {
    unsigned hi = 0;
    for (auto* obj : *dir->GetListOfKeys()) {
      TKey *key = static_cast<TKey*>(obj);
      TClass *key_class = TClass::GetClass(key->GetClassName());
      if (key_class->InheritsFrom(TH1::Class())) { // is TH1
        TH1 *h = static_cast<TH1*>(key->ReadObj());
        const int nbins = h->GetNbinsX()+2;
        auto it = std::find_if(hists.begin(),hists.end(),
          [h](const auto& pair){
            return !strcmp(pair.first->GetName(),h->GetName());
          });

        if (it==hists.end()) { // reserve memory
          hists.emplace_back(h,decltype(it->second)(nbins));
          it = --hists.end();
          for (auto& v : it->second) v.resize(dirs.size());
        }
        auto& bin = it->second;
        for (int i=0; i<nbins; ++i) {
          bin[i][hi] = h->GetBinContent(i);
        }
      } // end is TH1
      ++hi;
    }
  }
  unsigned count = 0;
  for (const auto& hh : hists) { // check hists count
    if (!count) count = hh.second.size();
    else if (hh.second.size() != count) throw std::runtime_error(cat(
      "different count of histogram \'",hh.first,
      "\' (",hh.second.size()," instead of ",count,')'
    ));
  }
  return hists;
}

auto scale_unc(const hists_t& hists) {
  struct hvar { TH1F *up, *down; };
  std::vector<hvar> hs;
  hs.reserve(hists.size());

  for (const auto& hh : hists) {
    const TH1 *h1 = hh.first;
    TH1F *h2 = new TH1F(
      h1->GetName(), h1->GetTitle(),
      h1->GetNbinsX(), h1->GetXaxis()->GetXbins()->GetArray()
    );
    hs.push_back({h2,static_cast<TH1F*>(h2->Clone())});
    const auto& hb = hs.back();
    for (const auto& h : hh.second) {
      for (int i=h.size(); i; ) { --i;
        const auto minmax = std::minmax_element(h.begin(),h.end());
        (*hb.up  )[i] = *minmax.second;
        (*hb.down)[i] = *minmax.first;
      }
    }
  }
  return hs;
}

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "usage: " << argv[0] << " in.root out.root" << endl;
    return 1;
  }
  TH1::AddDirectory(false);

  TFile fin(argv[1]);
  if (fin.IsZombie()) return 1;
  TFile fout(argv[2],"recreate");
  if (fout.IsZombie()) return 1;

  dirs_t scales, pdfs;

  for (auto* obj : *fin.GetListOfKeys()) {
    TKey *key = static_cast<TKey*>(obj);
    TClass *key_class = TClass::GetClass(key->GetClassName());
    if (key_class->InheritsFrom(TDirectory::Class())) {
      const char *name = key->GetName();
      const char *scale_end = strstr(name,"_PDF");
      if (!scale_end) {
        cerr << "no PDF name in directory name" << endl;
        return 1;
      }
      const char *pdf_begin = scale_end+3;
      const char *pdf_end   = strchr(pdf_begin,'_');
      std::string scale_name(name,scale_end-name);
      std::string pdf_name(pdf_begin,pdf_end-pdf_begin);

      TDirectory *dir = static_cast<TDirectory*>(key->ReadObj());
      if (scales.size()==0) {
        scales.emplace_back(dir,std::move(scale_name));
        pdfs.emplace_back(dir,std::move(pdf_name));
      } else if (scale_name!=scales.front().name) {
        if (std::find_if(pdfs.begin(),pdfs.end(),
            [&](const auto& p){ return p.name == pdf_name; })!=pdfs.end()) {
          cerr << "repeated PDF name in " << name << endl;
          return 1;
        }
        scales.emplace_back(dir,std::move(scale_name));
      } else if (pdf_name!=pdfs.front().name) {
        if (std::find_if(scales.begin(),scales.end(),
            [&](const auto& p){ return p.name == scale_name; })!=scales.end()) {
          cerr << "repeated scale name in " << name << endl;
          return 1;
        }
        pdfs.emplace_back(dir,std::move(pdf_name));
      } else {
        cerr << "repeated scale and PDF name in " << name << endl;
        return 1;
      }
    }
  }

  const auto _scale_unc = scale_unc(get_values(scales));
  const auto   _pdf_unc =   pdf_unc(get_values(pdfs));
  
}
