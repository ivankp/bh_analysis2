#ifndef NLO_MULTIBIN_HH
#define NLO_MULTIBIN_HH

#include <vector>
#include "multibin.hh"

struct nlo_multibin: public multibin {
  struct bin {
    int id = 0;
    double wtmp = 0, w = 0, w2 = 0;
  };
  size_t n = 0;
  std::vector<bin> bins;
  nlo_multibin(): bins(weights.size()) { }

  inline nlo_multibin& operator++() noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const double weight = weights[i];
      bin& b = bins[i];
      if (b.id == current_id) b.wtmp += weight;
      else {
        b.id = current_id;
        b.w2 += b.wtmp*b.wtmp;
        b.wtmp = weight;
      }
      b.w += weight;
    }
    ++n;
    return *this;
  }
  inline nlo_multibin& operator+=(const nlo_multibin& rhs) noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const bin& br = rhs.bins[i];
      bin& bl = bins[i];
      bl.wtmp += br.wtmp;
      bl.w += br.w;
      bl.w2 += br.w2;
    }
    n += rhs.n;
    return *this;
  }
};

namespace ivanp { namespace root {
template <> class bin_converter<nlo_multibin> {
  inline const auto& get(const nlo_multibin& bin) const noexcept {
    return bin.bins[nlo_multibin::wi];
  }
public:
  inline auto weight(const nlo_multibin& b) const noexcept {
    return get(b).w;
  }
  inline auto sumw2(const nlo_multibin& b) const noexcept {
    auto _b = get(b);
    return _b.w2 + _b.wtmp*_b.wtmp;
  }
  inline auto num(const nlo_multibin& b) const noexcept {
    return b.n;
  }
};
}}

#endif
