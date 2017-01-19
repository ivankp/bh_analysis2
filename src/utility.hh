#ifndef IVANP_UTILITY_HH
#define IVANP_UTILITY_HH

#include <string>

namespace ivanp {

template <typename T>
struct named {
  using type = T;
  type *p;
  std::string name;

  named(): p(nullptr), name() { }
  named(const named& n) = default;
  named(named&& n) = default;
  template <typename P, typename N>
  named(P&& ptr, N&& name)
  : p(std::forward<P>(ptr)), name(std::forward<N>(name)) { }

  inline type& operator*() const noexcept { return *p; }
  inline type* operator->() const noexcept { return p; }
};

#ifdef _GLIBCXX_VECTOR
template <typename T>
auto reserve(size_t n) {
  std::vector<T> _v;
  _v.reserve(n);
  return _v;
}
#endif

}

#endif
