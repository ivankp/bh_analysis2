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
#include <TAxis.h>

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
using hists_t = std::vector<std::pair<TH1*,std::vector<std::vector<double>>>>;

class {
private:
  bool on = false;
  double scale = 1.;
public:
  inline void set(double s) noexcept { on = true; scale = s; }
  inline double operator()(double x) const noexcept { return on ? x*scale : x; }
  inline double operator*() const noexcept { return scale; }
} scale;

auto get_values(const dirs_t& dirs) {
  hists_t hists;
  for (const auto& dir : dirs) {
    auto& keys = *dir->GetListOfKeys();
    for (auto* obj : keys) {
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
          for (auto& v : it->second) v.reserve(keys.GetSize());
        }
        if (h->Integral()==0.) h = it->first; // hack to remove black histograms
        for (int i=0; i<nbins; ++i) {
          it->second[i].push_back(scale(h->GetBinContent(i)));
        }
      } // end is TH1
    }
  }
  unsigned count = 0;
  for (const auto& h : hists) { // check hists count
    for (unsigned i=0, n=h.second.size(); i<n; ++i) {
      const auto& bin = h.second[i];
      if (!count) count = bin.size();
      else if (bin.size() != count) throw std::runtime_error(cat(
        "different count in bin ",i,
        " of histogram \'",h.first->GetName(),
        "\' (",bin.size()," instead of ",count,')'
      ));
    }
  }
  return hists;
}

auto mkhist(const TH1* h) {
  const auto* ax = h->GetXaxis();
  const auto* edges = ax->GetXbins()->GetArray();
  return edges
    ? new TH1F( h->GetName(), h->GetTitle(), h->GetNbinsX(), edges )
    : new TH1F( h->GetName(), h->GetTitle(), h->GetNbinsX(),
                ax->GetXmin(), ax->GetXmax() );
}
// template <typename T>
// inline auto clone(T* h) { return static_cast<T*>(h->Clone()); }

template <bool=true> struct envelope {
  TH1 *central;
  TH1F *up, *down;
  envelope(TH1* h): central(h), up(mkhist(h)), down(mkhist(h)) { }
};
template <> struct envelope<false> {
  TH1F *up, *down;
  envelope(TH1* h): up(mkhist(h)), down(mkhist(h)) { }
};

template <bool with_central=true>
auto mk_envelope(const hists_t& hists) {
  std::vector<envelope<with_central>> hs;
  hs.reserve(hists.size());
  for (const auto& hh : hists) {
    hs.emplace_back(hh.first);
    const auto& hb = hs.back();
    for (int i=hh.second.size(); i; ) { --i;
      const auto minmax = std::minmax_element(
        hh.second[i].begin(),hh.second[i].end());
      (*hb.up  )[i] = *minmax.second - hh.second[i][0];
      (*hb.down)[i] = hh.second[i][0] - *minmax.first;
    }
  }
  return hs;
}

auto mk_pdf_unc(const hists_t& hists, const LHAPDF::PDFSet& pdf) {
  struct hvar { TH1F *errplus, *errminus, *errsymm; };
  std::vector<hvar> hs;
  hs.reserve(hists.size());
  for (const auto& hh : hists) {
    const TH1 *h = hh.first;
    hs.push_back({mkhist(h),mkhist(h),mkhist(h)});
    const auto& hb = hs.back();
    for (int i=hh.second.size(); i; ) { --i;
      const auto unc = pdf.uncertainty(hh.second[i]);
      (*hb.errplus )[i] = /*hh.second[i][0] +*/ unc.errplus;
      (*hb.errminus)[i] = /*hh.second[i][0] -*/ unc.errminus;
      (*hb.errsymm )[i] = unc.errsymm;
    }
  }
  return hs;
}

int main(int argc, char* argv[]) {
  if (argc!=3 && argc!=4) {
    cout << "usage: " << argv[0] << " in.root out.root scale" << endl;
    return 1;
  }
  if (argc==4) { scale.set(atof(argv[3])); }

  TH1::AddDirectory(false);

  TFile fin(argv[1]);
  if (fin.IsZombie()) return 1;

  dirs_t scales, pdfs;
  std::string central_scale_name, central_pdf_name;

  for (auto* obj : *fin.GetListOfKeys()) {
    TKey *key = static_cast<TKey*>(obj);
    TClass *key_class = TClass::GetClass(key->GetClassName());
    if (key_class->InheritsFrom(TDirectory::Class())) {
      const char *name = key->GetName();
      cout << name << endl;
      const char *scale_end = strstr(name,"_PDF");
      if (!scale_end) {
        cerr << "no PDF name in directory name" << endl;
        return 1;
      }
      const char *pdf_begin = scale_end+4;
      const char *pdf_end   = strchr(pdf_begin,'_');
      std::string scale_name(name,scale_end-name);
      std::string pdf_name(pdf_begin,pdf_end-pdf_begin);

      TDirectory *dir = static_cast<TDirectory*>(key->ReadObj());
      if (scales.size()==0) {
        central_scale_name = scale_name;
        central_pdf_name = pdf_name;
        scales.emplace_back(dir,std::move(scale_name));
        pdfs.emplace_back(dir,std::move(pdf_name));
      } else if (scale_name!=central_scale_name && pdf_name==central_pdf_name) {
        if (std::find_if(scales.begin(),scales.end(),
            [&](const auto& p){ return p.name == scale_name; })!=scales.end()) {
          cerr << "repeated scale name in " << name << endl;
          return 1;
        }
        scales.emplace_back(dir,std::move(scale_name));
      } else if (scale_name==central_scale_name && pdf_name!=central_pdf_name) {
        if (std::find_if(pdfs.begin(),pdfs.end(),
            [&](const auto& p){ return p.name == pdf_name; })!=pdfs.end()) {
          cerr << "repeated PDF name in " << name << endl;
          return 1;
        }
        pdfs.emplace_back(dir,std::move(pdf_name));
      } else if (scale_name!=central_scale_name && pdf_name!=central_pdf_name) {
        cerr << "unexpected scale and PDF name in " << name << endl;
        return 1;
      } else {
        cerr << "repeated scale and PDF name in " << name << endl;
        return 1;
      }
    }
  }

  TFile fout(argv[2],"recreate","histograms with uncertainties",109);
  if (fout.IsZombie()) return 1;

  const auto pdf = LHAPDF::lookupPDF(std::stoi(central_pdf_name));
  if (pdf.second==-1) throw std::runtime_error(
    "cannot find PDF for "+central_pdf_name);
  cout << "\nPDF set: " << pdf.first << endl;

  const auto scale_unc = mk_envelope(get_values(scales));
  const auto   pdf_unc = mk_pdf_unc(get_values(pdfs), LHAPDF::PDFSet(pdf.first));
  // const auto   pdf_unc =   mk_envelope<false>(get_values(pdfs));

  auto *dir = fout.mkdir("central");
  dir->cd();
  for (const auto& h : scale_unc) {
#if ROOT_VERSION_CODE <= 394760
    // hack to remove bin labels
    memset(reinterpret_cast<char*>(h.central->GetXaxis())
           + sizeof(TAxis) - sizeof(void*), 0, sizeof(void*));
#endif
    h.central->Scale(*scale);
    h.central->SetDirectory(dir);
  }

#define VAR(TYPE,NAME) \
  dir = fout.mkdir(#TYPE"_"#NAME); \
  dir->cd(); \
  for (const auto& h : TYPE##_unc) h.NAME->SetDirectory(dir);

  VAR(scale,up)
  VAR(scale,down)
  // VAR(pdf,up)
  // VAR(pdf,down)
  VAR(pdf,errplus)
  VAR(pdf,errminus)
  VAR(pdf,errsymm)

  fout.Write();
}
