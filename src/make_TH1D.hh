#ifndef IVANP_BINNER_MAKE_TH1D
#define IVANP_BINNER_MAKE_TH1D

template <typename A, typename B>
TH1D* make_TH1D(const char* name, const A& axis, B begin, B end) {
  TH1D* h = new TH1D(name,"",axis.nbins(),axis.min(),axis.max());
  h->Sumw2();
  TArrayD& sumw2 = *h->GetSumw2();
  size_t n_total = 0, i = 0;
  for (auto bin=begin; bin!=end; ++bin, ++i) {
    h->SetBinContent(i,bin->w);
    sumw2[i] = bin->w2;
    n_total += bin->n;
  }
  h->SetEntries(n_total);
  return h;
}

template <typename A>
TH1D* root_hist(const hist_t<A>& h, const std::string& name) {
  return make_TH1D(name.c_str(),h.axis(),h.bins().begin(),h.bins().end());
}

template <typename A1, typename A2>
void root_hist(const hist_t<A1,A2>& h, const std::string& name,
  const std::string var2
) {
  // TODO: compiler bug?
  // const auto& a1 = h.axis<0>();
  // const auto& a2 = h.axis<1>();

  const auto& a1 = std::get<0>(h.axes());
  const auto& a2 = std::get<1>(h.axes());
  const auto nbins1 = a1.nbins()+2;
  const auto nbins2 = a2.nbins()+1;
  for (const auto& b : h.bins()) cout << ' ' << b.w;
  cout << endl;
  test(h.bins().size())
  test(nbins1)
  test(nbins2)
  auto it = h.bins().begin() + nbins1;
  for (unsigned i=1; i<nbins2; ++i) {
    make_TH1D(cat(name,'_',var2,"_[",a2.lower(i),',',a2.upper(i),')').c_str(),
      a1,it,it+nbins1);
    it += nbins1;
  }
}

#endif
