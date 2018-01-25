#ifndef BIN_DEFS_HH
#define BIN_DEFS_HH

#include <vector>

template <typename T>
struct arrow_self {
  inline T* operator->() noexcept {
    return static_cast<T*>(this);
  }
  inline const T* operator->() const noexcept {
    return static_cast<const T*>(this);
  }
};

struct lo_bin: arrow_self<lo_bin> {
  double w = 0, w2 = 0;
  long unsigned n = 0;
  inline void operator()() noexcept {
    ++w; ++w2; ++n;
  }
  inline void operator()(double weight) noexcept {
    w  += weight;
    w2 += weight*weight;
    ++n;
  }
  inline lo_bin& operator+=(const lo_bin& b) noexcept {
    w  += b.w;
    w2 += b.w2;
    n  += b.n;
    return *this;
  }
  inline double err2() const noexcept { return w2; }
};

struct nlo_bin: arrow_self<nlo_bin> {
  int id = 0; // event id
  static int current_id;
  double w = 0, wtmp = 0, w2 = 0;
  long unsigned n = 0;
  inline void operator()() noexcept {
    if (id == current_id) ++wtmp;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      ++wtmp;
    }
    ++w;
    ++n;
  }
  inline void operator()(double weight) noexcept {
    if (id == current_id) wtmp += weight;
    else {
      id = current_id;
      w2 += wtmp*wtmp;
      wtmp = weight;
    }
    w += weight;
    ++n;
  }
  inline nlo_bin& operator+=(const nlo_bin& b) noexcept {
    wtmp += b.wtmp;
    w  += b.w;
    w2 += b.w2;
    n  += b.n;
    return *this;
  }
  inline double err2() const noexcept { return w2 + wtmp*wtmp; }
};
int nlo_bin::current_id;

struct profile_bin: arrow_self<profile_bin> {
  double w = 0, // total weight
         m = 0, // mean
         s = 0; // variance * (n-1)
  long unsigned n = 0;
  inline void operator()(double x) noexcept {
    ++w;
    if (n == 0) {
      m = x;
    } else {
      const double d = x-m;
      m += d/w;
      s += d*(x-m);
    }
    ++n;
  }
  inline void operator()(double weight, double x) noexcept {
    w += weight;
    if (n == 0) {
      m = x;
    } else {
      const double wd = weight*(x-m);
      m += wd/w;
      s += wd*(x-m);
    }
    ++n;
  }
  inline double var() const noexcept { return (n > 1) ? n*s/((n-1)*w) : 0.; }
};

struct multiweight_bin_base {
  static std::vector<double> weights;
  static unsigned wi;
};
std::vector<double> multiweight_bin_base::weights;
unsigned multiweight_bin_base::wi;

template <typename Bin>
struct multiweight_bin: multiweight_bin_base {
  std::vector<Bin> bins;
  multiweight_bin(): bins(weights.size()) { }

  template <typename... T>
  inline multiweight_bin& operator()(T... x) noexcept {
    for (unsigned i=weights.size(); i; ) {
      --i; bins[i](weights[i],x...);
    }
    return *this;
  }
  inline multiweight_bin& operator+=(const multiweight_bin& rhs) noexcept {
    for (unsigned i=weights.size(); i; ) {
      --i; bins[i] += rhs.bins[i];
    }
    return *this;
  }
  inline Bin& operator->() noexcept { return bins[wi]; }
  inline const Bin& operator->() const noexcept { return bins[wi]; }
};

namespace ivanp { namespace root {
template <typename Bin> struct bin_converter<multiweight_bin<Bin>> {
  using type = multiweight_bin<Bin>;
  inline auto val (const type& b) const noexcept { return b->w; }
  inline auto err2(const type& b) const noexcept { return b->err2(); }
  inline auto num (const type& b) const noexcept { return b->n; }
};
template <> struct bin_converter<profile_bin> {
  using type = profile_bin;
  inline auto val (const type& b) const noexcept { return b->m; }
  inline auto err2(const type& b) const noexcept { return b->var(); }
  inline auto num (const type& b) const noexcept { return b->n; }
};
template <> struct bin_converter<multiweight_bin<profile_bin>> {
  using type = multiweight_bin<profile_bin>;
  inline auto val (const type& b) const noexcept { return b->m; }
  inline auto err2(const type& b) const noexcept { return b->var(); }
  inline auto num (const type& b) const noexcept { return b->n; }
};
}}

#endif
