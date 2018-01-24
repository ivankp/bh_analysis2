#ifndef NLO_MULTIBIN_HH
#define NLO_MULTIBIN_HH

#include <vector>
#include "multibin.hh"

struct nlo_bin {
  int id = 0;
  static int current_id;
  double wtmp = 0, w = 0, w2 = 0;
  inline nlo_bin& operator+=(double weight) noexcept {
    if (id == current_id) wtmp += weight;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      wtmp = weight;
    }
    w += weight;
    return *this;
  }
  inline nlo_bin& operator+=(const nlo_bin& b) noexcept {
    wtmp += b.wtmp;
    w += b.w;
    w2 += b.w2;
    return *this;
  }
  inline double sumw2() const noexcept { return w2 + wtmp*wtmp; }
  inline nlo_bin* operator->() noexcept { return this; }
  inline const nlo_bin* operator->() const noexcept { return this; }
};
int nlo_bin::current_id;

template <typename Bin = nlo_bin>
struct nlo_multibin: public multibin {
  size_t n = 0;
  std::vector<Bin> bins;
  nlo_multibin(): bins(weights.size()) { }

  inline nlo_multibin& operator++() noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      bins[i] += weights[i];
    }
    ++n;
    return *this;
  }
  inline nlo_multibin& operator+=(const nlo_multibin& rhs) noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      bins[i] += rhs.bins[i];
    }
    n += rhs.n;
    return *this;
  }
  inline Bin& operator->() noexcept { return bins[wi]; }
  inline const Bin& operator->() const noexcept { return bins[wi]; }
};

namespace ivanp { namespace root {
template <typename Bin> struct bin_converter<nlo_multibin<Bin>> {
  using type = nlo_multibin<Bin>;
  inline auto weight(const type& b) const noexcept { return b->w; }
  inline auto sumw2 (const type& b) const noexcept { return b->sumw2(); }
  inline auto num   (const type& b) const noexcept { return b.n; }
};
}}

#endif
