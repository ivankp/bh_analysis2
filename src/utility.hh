#ifndef IVANP_UTILITY_HH
#define IVANP_UTILITY_HH

#include <string>

namespace ivanp {

template <typename T>
struct named_ptr {
  using type = T;
  type *p;
  std::string name;

  named_ptr(): p(nullptr), name() { }
  named_ptr(const named_ptr& n) = default;
  named_ptr(named_ptr&& n) = default;
  template <typename P, typename N>
  named_ptr(P&& ptr, N&& name)
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
