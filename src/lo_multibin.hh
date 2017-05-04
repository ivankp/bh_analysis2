#ifndef LO_MULTIBIN_HH
#define LO_MULTIBIN_HH

#include "multibin.hh"

struct lo_multibin: public multibin {
  struct bin {
    double w = 0., w2 = 0.;
  };
  size_t n = 0;
  std::vector<bin> bins;
  lo_multibin(): bins(weights.size()) { }

  inline lo_multibin& operator++() noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const double weight = weights[i];
      bin& b = bins[i];
      b.w  += weight;
      b.w2 += weight*weight;
    }
    ++n;
    return *this;
  }
  inline lo_multibin& operator+=(const lo_multibin& rhs) noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const bin& br = rhs.bins[i];
      bin& bl = bins[i];
      bl.w  += br.w;
      bl.w2 += br.w2;
    }
    n += rhs.n;
    return *this;
  }
};

namespace ivanp { namespace root {
template <> class bin_converter<lo_multibin> {
  inline const auto& get(const lo_multibin& bin) const noexcept {
    return bin.bins[lo_multibin::wi];
  }
public:
  inline auto weight(const lo_multibin& b) const noexcept { return get(b).w;  }
  inline auto sumw2 (const lo_multibin& b) const noexcept { return get(b).w2; }
  inline auto num   (const lo_multibin& b) const noexcept { return b.n; }
};
}}

#endif
