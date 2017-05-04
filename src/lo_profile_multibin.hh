#ifndef LO_PROFILE_MULTIBIN_HH
#define LO_PROFILE_MULTIBIN_HH

#include "multibin.hh"

struct lo_profile_multibin: public multibin {
  struct bin {
    double w = 0., m1, m2, s1, s2;
  };
  size_t n = 0;
  std::vector<bin> bins;
  lo_profile_multibin(): bins(weights.size()) { }

  inline lo_profile_multibin& operator()(double x) noexcept {
    for (unsigned i=weights.size(); i!=0; ) {
      --i;
      const double weight = weights[i];
      bin& b = bins[i];
      b.w += weight;
      const double xw = x*weight;
      if (n == 0) {
        b.m1 = b.m2 = xw;
        b.s1 = 0.;
      } else {
        const double xw_m1w = xw - weight*b.m1;
        b.m2 = b.m1 + xw_m1w/b.w;
        b.s2 = b.s1 + xw_m1w*(xw - weight*b.m2);
        // set up for next iteration
        b.m1 = b.m2; 
        b.s1 = b.s2;
      }
    }
    ++n;
    return *this;
  }
};

namespace ivanp { namespace root {
template <> class bin_converter<lo_profile_multibin> {
  inline const auto& get(const lo_profile_multibin& bin) const noexcept {
    return bin.bins[lo_profile_multibin::wi];
  }
public:
  inline auto weight(const lo_profile_multibin& b) const noexcept {
    return b.n ? get(b).m2 : 0.;
  }
  inline auto sumw2 (const lo_profile_multibin& b) const noexcept {
    auto _b = get(b);
    return (_b.w > 1) ? _b.s2/(_b.w - 1) : 0.;
  }
  inline auto num(const lo_profile_multibin& b) const noexcept { return b.n; }
};
}}

#endif
